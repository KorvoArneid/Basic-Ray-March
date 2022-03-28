#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cassert>
#include <random>
#include <chrono>
#include <utility>
#include <memory>

using namespace std;

using flt=double; //precision and speed balance
const flt inf=1e10;
const flt eps=1e-6;
template<typename T>
using ptr=shared_ptr<T>;

flt randflt(){ //returns float number between 0 and 1 inclusive. Adjust r variable for precision;
    int l=0, r=1000000;
    int val;
    static mt19937_64 gen(chrono::steady_clock::now().time_since_epoch().count()); val=uniform_int_distribution<int>(l, r)(gen);
    return val*1. / r;
}

struct point{
    flt x, y, z;
    point(flt x=0, flt y=0, flt z=0): x(x), y(y), z(z) {}
    friend point operator+(point a, point b){
        return point(a.x+b.x, a.y+b.y, a.z+b.z);
    }
    friend point operator-(point a, point b){
        return point(a.x-b.x, a.y-b.y, a.z-b.z);
    }
    friend point operator+(point a, flt b){
        return point(a.x+b, a.y+b, a.z+b);
    }
    friend point operator-(point a, flt b){
        return point(a.x-b, a.y-b, a.z-b);
    }
    friend point operator*(point a, flt b){
        return point(a.x*b, a.y*b, a.z*b);
    }
    friend point operator/(point a, flt b){
        b=1./b;
        return a*b;
    }
    void operator+=(point b){
        *this=*this+b;
    }
    void operator+=(flt b){
        *this=*this+b;
    }
    void operator-=(point b){
        *this=*this-b;
    }
    void operator-=(flt b){
        *this=*this-b;
    }
    void operator*=(flt b){
        *this=*this*b;
    }
    void operator/=(flt b){
        *this=*this/b;
    }
    flt dot(point b) const{
        return x*b.x + y*b.y + z*b.z;
    }
    point cross(point b) const{
        flt px=y*b.z-b.y*z;
        flt py=x*b.z-b.x*z;
        flt pz=x*b.y-b.x*y;
        auto ret=point(px, -py, pz);
        ret.normalize();
        return ret;
    }
    flt invLen() const{  //quake inverse sqrt trick
        flt number=x*x+y*y+z*z;
        long i;
	    float x2, y;
	    const float threehalfs = 1.5F;
        x2 = number * 0.5F;
	    y  = number;
	    i  = * ( long * ) &y;                       
	    i  = 0x5f3759df - ( i >> 1 );              
	    y  = * ( float * ) &i;
	    y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 3rd iteration??!!
	    return y;
    }
    flt len() const{
        return sqrt(x*x + y*y + z*z);
    }
    flt flen() const{ //fast length using quake inverse sqrt technique. Might have precision issues
       return 1./invLen();
    }
    void normalize(){
        flt inv=invLen();
        if(inv>inf) {
            cerr<<"Inverse is greater than infinite in normalize\n";
        }
        *this = point(x*inv, y*inv, z*inv);
    }
    friend ostream& operator<<(ostream& out, point a){
        out<<"("<<a.x<<", "<<a.y<<", "<<a.z<<")";
        return out;
    }
};

struct ray{
    point o, dir;
    ray(point o=point(), point d=point()): o(o), dir(d) { d.normalize(); };
};

struct object{
    vector<ptr<object>> pars;
    int mode=1000;
    private:
    object(const object& copy){} //can't copy
    public:
    object() {}
    object(const ptr<object>& a){
        pars.push_back(a);
        mode=0;
    }
    void add(const ptr<object>& b){
        mode=1;
        if(pars.size()==1) pars.push_back(b);
        else pars[1]=b;
    }
    void subtract(ptr<object> b){
        mode=-1;
        if(pars.size()==1) pars.push_back(b);
        else pars[1]=b;
    }
    void remove() {//removes and subtracted and added objects
        mode=0;
    }
    virtual flt dist(point p) const{
        if(mode==0){
            return pars[0]->dist(p);
        }
        else if(mode==1){
            return min(pars[0]->dist(p), pars[1]->dist(p));
        }
        else if(mode==-1){
            return max(pars[0]->dist(p), -pars[1]->dist(p));
        }
        assert(string("Wrong mode in object class").length()==0);
        return 0;
    }
    virtual bool inside(point p) const{
        if(mode==0){
            return pars[0]->inside(p);
        }
        else if(mode==1){
            return (pars[0]->inside(p) | pars[1]->inside(p));
        }
        else if(mode==-1){
            bool in1=pars[0]->inside(p);
            if(!in1) return 0;
            if(pars[1]->inside(p)) return 0;
            else return 1;
        }
        assert(string("Wrong mode in object class").length()==0);
        return 0;
    }
    virtual point surfaceSample() const{
        int count=0;
        int maxlimit = 10; //how many random samples before
        if(mode==0) return pars[0]->surfaceSample();
        point ans;
        while(count++ < maxlimit){
            point sur1=pars[0]->surfaceSample();
            point sur2=pars[1]->surfaceSample();
            if(mode==1){
                if(!pars[1]->inside(sur1)){
                    ans=sur1;
                    break;
                }
                else if(!pars[0]->inside(sur2)){
                    ans=sur2;
                    break;
                }
            }
            else if(mode==-1){ //mode 2 / subtractive
                if(!pars[1]->inside(sur1)) {
                    ans=sur1;
                    break;
                }
                else if(pars[0]->inside(sur2)){
                    ans=sur2;
                    break;
                }
            }
        }
        assert(mode!=1000); //1000 mode means empty object 
        assert(count<=maxlimit);
        return ans;
    }
    virtual bool onSurface(point p) const{
        if(mode==0) return pars[0]->onSurface(p);
        else if(mode==1){
            if(pars[0]->onSurface(p) && !pars[1]->inside(p)) return 1;
            if(pars[1]->onSurface(p) && !pars[0]->inside(p)) return 1;
            return 0;
        }
        else if(mode==-1){
            if(pars[0]->onSurface(p) && !pars[1]->inside(p)) return 1;
            if(pars[0]->inside(p) && pars[1]->onSurface(p)) return 1;
            return 0;
        }
        assert(string("Wrong mode in object class").length()==0);
        return 0;
    }
    virtual point normal(point p) const{
        if(mode==0) return pars[0]->normal(p);
        else if(mode==1){
            if(pars[0]->onSurface(p)){
                assert(pars[1]->onSurface(p) || !pars[1]->inside(p));
                return pars[0]->normal(p);
            }
            else if(pars[1]->onSurface(p)){
                assert(pars[0]->onSurface(p) || !pars[0]->inside(p));
                return pars[1]->normal(p);
            }
            else assert(string("Not on surface of object").length()==0);
        }
        else if(mode==-1){
            if(pars[0]->onSurface(p) && !pars[1]->inside(p)) return pars[0]->normal(p);
            else if(pars[1]->onSurface(p) && pars[0]->inside(p)) return pars[1]->normal(p) * -1; //reversing the normal
            else assert(string("Not on surface of object").length()==0);
        }
        assert(string("Wrong mode in object class").length()==0);
        return 0;
    }
};

/*object subclasses follow pattern:
class sub{
   virtual flt dist(point p) const;
   virtual bool inside(point p) const;
   virtual point surfaceSample() const;
   virtual bool onSurface(point p) const;
   virtual point normal(point p) const;
}*/

struct sphere : public object{
    point o; flt rad;
    sphere(point o, flt rad): o(o), rad(rad) {}
    virtual flt dist(point p) const{
        return (o-p).flen()-rad;
    }
    virtual bool inside(point p) const{
        return (dist(p)<eps); //eps or 0?
    }
    virtual point surfaceSample() const{
        point vec(2*randflt()-1, 2*randflt()-1, 2*randflt()-1);
        vec.normalize();
        return o + (vec*rad);
    }
    virtual bool onSurface(point p) const{
        return (abs(dist(p))<eps);
    }
    virtual point normal(point p) const{
        assert(onSurface(p));
        p=(p-o);
        p.normalize();
        return p;
    }
};

struct lineSegment : public object {
    point a, b;
    flt rad;
    lineSegment(point a, point b, flt rad) : a(a), b(b), rad(rad) {}
    private:
    point closest_point(point p) const{
        flt h=min(1.0, max(0.0, (p-a).dot(b-a) / (b-a).dot(b-a)));
        return a+(b-a)*h;
    }
    public:
    virtual flt dist(point p) const{
        return (p-closest_point(p)).len()-rad;
    }
    virtual bool inside(point p) const{
        return (dist(p)<eps);
    }
    virtual point surfaceSample() const{
        assert(1==0); //AHAHAH
    }
    virtual bool onSurface(point p) const{
        return (abs(dist(p))<eps);
    }
    virtual point normal(point p) const{
        assert(onSurface(p));
        p=p-closest_point(p);
        p.normalize();
        return p;
    }
};

struct camera{
    point pos;
    flt rot;
    double aspect;
    double width;
    double height;
    int row, col;
    double zoom=1;
    camera(point pos, double rot, int ro, int co, flt w=1, double z=1): pos(pos), rot(rot), width(w), row(ro), col(co), zoom(z) {
        aspect= row*1./col;
        height=width / aspect;
    }
    ray give_ray(int u, int v){ //u is vertical, v is horizontal
        point center=pos;
        center.z+=zoom;
        point ll(center.x-width/2, center.y-height/2, center.z);
        point ur(center.x+width/2, center.y+height/2, center.z);
        ray ret(pos, point(1, 1, 1));
        ret.dir=point(ll.x+(v*1./col)*(ur.x-ll.x), ll.y+(u*1./row)*(ur.y-ll.y), center.z);
        ret.dir-=pos;
        ret.dir.normalize();
        return ret;
    }
 //later :p
};

struct scene{
    vector<ptr<object>> objs;
    int count=0;
    scene(){}
    void add(const ptr<object>& p){
        objs.push_back(p);
        count++;
    }
    void add(const scene& s){
        for(auto it:s.objs) add(it);
    }
    friend scene operator+(const scene& a, const scene& b){
        scene tmp;
        for(auto it:a.objs) tmp.add(it);
        for(auto it:b.objs) tmp.add(it);
        return tmp;
    }
    pair<flt, ptr<object>> march(point p) const{
        flt ansd=inf;
        ptr<object> ansp;
        for(auto it : objs) {  //finding the march distance
            flt d=it->dist(p);
            if(ansd>d) {
                ansp=it;
                ansd=d;
            }
        }
        return make_pair(ansd, ansp);
    }  
};

struct engine{

    void adjustCol(point& col){
        flt maxBrightness=0.01;
        col.x=(col.x>maxBrightness)? 1 : ((col.x<0)? 0 : col.x/maxBrightness);
        col.y=(col.y>maxBrightness)? 1 : ((col.y<0)? 0 : col.y/maxBrightness);
        col.z=(col.z>maxBrightness)? 1 : ((col.z<0)? 0 : col.z/maxBrightness);
    }

    point collided(point p, const ptr<object>& it, const scene& lights, const scene& world){
        flt res=0;
        int count=0;
        point nor=it->normal(p);
        for(auto light:lights.objs){
            if(light==it){ //point is on a light sourrce itself
                res=1;
                count=1;
                break;
            }
            int samples=0;
            int maxSamples=10;
            while(samples++<maxSamples){
                point random=light->surfaceSample();
                ray r(p, random-p);
                ptr<object> tmp(new object);
                bool ok=march(world, r, tmp);
                flt dist=(random-p).flen();
                if(!ok || (r.o-p).flen() > dist){
                    flt dot=nor.dot((random-p)/dist);
                    res+=1./(dist*dist)*((dot<0)? 0 : dot*dot);
                }
            }
            count+=maxSamples;
        }
        point ret(res/count, res/count, res/count);
        adjustCol(ret);
        return ret;
    }

    bool march(const scene& world, ray& r, ptr<object>& which){
        point cur;
        bool outside=1;
        while((cur-point(0, 0, 0)).flen()<inf){
            auto mar=world.march(cur);
            assert(mar.first>=-eps); //doesn't enter objects
            if(mar.first<eps){
                which=mar.second;
                r.o=cur;
                outside=0;
                break;
            }
            cur+=r.dir*mar.first;
        }
        return !outside;
    }

    engine(ofstream& file, int width, int height, const scene& world, const scene& lights){ //no ray bounces
        camera cam(point(), 0, height, width);
        scene allobjs=world+lights;
        int Pixels=width*height;
        for(int u=height-1; u>=0; u--){
            for(int v=0; v<width; v++){
                ray r=cam.give_ray(u, v);
                ptr<object> tmp(new object); //temporary pointer to contain object with which collision occurs(if)
                if(!march(allobjs, r, tmp)) file<<"0 0 0\n";
                else {
                    point col=collided(r.o, tmp, lights, world);
                    file<<int(col.x*255)<<" "<<int(col.y*255)<<" "<<int(col.z*255)<<"\n";
                }
                //cout<<"\rProgress... "<<((height-1-u)*width+v+1)*100 / Pixels<<"%";  //uncomment if you want updates for each pixel
            }
            cout<<"\rProgress... "<<((height-u)*width)*100 / Pixels<<"%";
        }
    }
};

ptr<object> makeObj(const vector<pair<int, ptr<object>>>& objs){ //add objects ptrs in pair by adding 1 to add and -1 to subtract the object in the final object pointer 
    assert(objs.size()>0 && objs[0].first==1);
    ptr<object> ret=objs[0].second;
    for(int i=1; i<(int)objs.size(); i++){
        if(objs[i].first==1){
            ret=ptr<object>(new object(ret));
            ret->add(objs[i].second);
        }
        else if(objs[i].first==-1){
            ret=ptr<object>(new object(ret));
            ret->subtract(objs[i].second);
        }
        else {
            assert(string("Wrong choice for object addition or subtraction!!").length()==0);
        }
    }
    return ret;
}

int main(){
    //freopen("Log.txt", "w", stderr);
    auto file = ofstream("image.ppm");
    int width=1000, height=1000;
    file<<"P3\n"<<width<<" "<<height<<"\n255\n";
    scene world;
    scene lights;
    //ptr<T> foo(new foo(<constructor params>));
    {
        //Cry emoji
        vector<pair<int, ptr<object>>> List;
        List.push_back(make_pair(1, ptr<object>(new object(ptr<object>(new sphere(point(0, 0, 5), 1.7)))))); //main face
        List.push_back(make_pair(1, ptr<object>(new object(ptr<object>(new sphere(point(-0.6, 0.5, 3.4), 0.3)))))); //left eyeball
        List.push_back(make_pair(1, ptr<object>(new object(ptr<object>(new sphere(point(0.6, 0.5, 3.4), 0.3)))))); //right eyeball
        List.push_back(make_pair(-1, ptr<object>(new object(ptr<object>(new sphere(point(-0.6, 0.5, 3.1), 0.1)))))); //left pupil
        List.push_back(make_pair(-1, ptr<object>(new object(ptr<object>(new sphere(point(0.6, 0.5, 3.1), 0.1)))))); //right pupil
        List.push_back(make_pair(-1, ptr<object>(new object(ptr<object>(new sphere(point(0, -0.5, 3.2), 0.3)))))); //mouth
        List.push_back(make_pair(1, ptr<object>(new object(ptr<object>(new lineSegment(point(-0.1, -0.55, 3.35), point(0.1, -0.55, 3.35), 0.1)))))); //tongue
        List.push_back(make_pair(-1, ptr<object>(new object(ptr<object>(new lineSegment(point(-0.6, -0.5, 3.4), point(-0.6, 0.2, 3.4), 0.2)))))); //left tear
        List.push_back(make_pair(-1, ptr<object>(new object(ptr<object>(new lineSegment(point(0.6, -0.5, 3.4), point(0.6, 0.2, 3.4), 0.2)))))); //right tear
        world.add(makeObj(List));
        lights.add(ptr<object>(new object(ptr<object>(new sphere(point(2, 2, 1), 1)))));
    }
    engine(file, width, height, world, lights);
    cout<<"\nDone!"<<endl;
}
