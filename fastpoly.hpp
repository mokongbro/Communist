#ifndef VEC_POLY_HPP
#define VEC_POLY_HPP
#include<vector>
#include<complex>
#include<compare>
#include<map>
#include<string>
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<type_traits>
#include<limits>
using i128=__int128;
using u128=unsigned __int128;
using u32=unsigned;
using u64=unsigned long long;
#define FPTP template<typename T>
using WO=std::weak_ordering;
FPTP using cx=std::complex<T>;
using Z=size_t;
constexpr u64 E18=1000000000000000000ull;
FPTP struct cxt:std::false_type{};
template<typename U>struct cxt<std::complex<U>>:std::true_type{};
FPTP using V=std::vector<T>;
using std::string;
using std::move;
using std::is_same_v;
FPTP constexpr bool fpt=std::is_floating_point_v<T>;
namespace fpf{
FPTP T pw10(int e){
	T r{1},b{10};
	if(e<0){e=-e;for(;e;e>>=1){if(e&1)r/=b;if(e>1)b*=b;}}
	else for(;e;e>>=1){if(e&1)r*=b;if(e>1)b*=b;}
	return r;
}
FPTP T cabs(const cx<T>&z){long double a{static_cast<long double>(z.real())},b{static_cast<long double>(z.imag())};return std::sqrt(a*a+b*b);}
FPTP T abs(const cx<T>&x){return cabs(x);}
FPTP T abs(const T&x){return x<0?-x:x;}
FPTP T trunc(const T&x){
	if(x!=x)return x;
	if constexpr(is_same_v<T,__float128>){
		if(x<0)return -trunc(-x);
		if(x>=1e38)return x;
		return static_cast<__float128>(static_cast<i128>(x));
	}else return std::trunc(x);
}
FPTP T mod(const T&a,const T&b){return a-b*trunc(a/b);}
FPTP bool eq(const T&x,const T&y,__float128 e){
	if(x==y)return true;
	if(!(x==x)||!(y==y))return false;
	T ax{abs(x)},ay{abs(y)};
	if((ax>0&&ax/2==ax)||(ay>0&&ay/2==ay))return false;
	T d{abs(x-y)},m{ax>ay?ax:ay};
	return d<=e||(m>1&&d<=e*m);
}
FPTP T num(const string&s){
	T v{};
	for(char c:s)if(c>='0'&&c<='9')v=v*10+(c-'0');
	return v;
}
}
FPTP WO operator<=>(const cx<T>&a,const cx<T>&b){
	T x{fpf::cabs(a)},y{fpf::cabs(b)};
	return x<y?WO::less:x>y?WO::greater:WO::equivalent;
}
FPTP WO operator<=>(const cx<T>&a,const typename cx<T>::value_type&b){return a<=>cx<T>{b};}
FPTP struct fpe{static constexpr __float128 v=1e-12;};
template<>struct fpe<__float128>{static constexpr __float128 v=1e-30;};
FPTP struct fpe<cx<T>>:fpe<T>{};
template<typename T=i128>
class Bint{
	public:
	static constexpr T base=E18;
	V<T>dt{};
	bool sig{};
	Bint():dt{0},sig{}{}
	Bint(const T&t):dt{},sig{t<0}{
		i128 x{t<0?-static_cast<i128>(t):static_cast<i128>(t)};
		put(*this,static_cast<u128>(x));
	}
	Bint(const string&s):dt{},sig{}{
		Z n{s.size()},d{};
		for(char c:s){if(c=='-'){sig=1;continue;}d++;}
		dt.assign((d+17)/18,{});
		Z k{};
		u64 w{},p{1};
		for(Z i{n};i--;){
			char c{s[i]};
			if(c=='-')continue;
			w+=(c-'0')*p;
			p*=10;
			if(p==E18){dt[k++]=w;w=0;p=1;}
		}
		if(p!=1)dt[k++]=w;
		reload();
	}
	void reload(){
		while(dt.size()>1&&!dt.back())dt.pop_back();
		if(dt.empty())dt.push_back({});
		if(dt.size()==1&&!dt[0])sig=0;
	}
	explicit operator bool()const{return dt.size()>1||dt[0];}
	static u64 div9(u128&x){
		u64 r{},lo=x,hi=x>>64;
		u64 w[4]{lo&4294967295ull,lo>>32,hi&4294967295ull,hi>>32},q[4]{};
		for(int i{3};i>=0;i--){u64 c{r*4294967296ull+w[i]};q[i]=c/1000000000ull;r=c-q[i]*1000000000ull;}
		x=(static_cast<u128>(q[3]*4294967296ull+q[2])<<64)|(q[1]*4294967296ull+q[0]);
		return r;
	}
	static u64 divbase(u128&x){
		u64 d0{div9(x)};
		u64 d1{div9(x)};
		return d0+d1*1000000000ull;
	}
	static void put(Bint&t,u128 x){
		t.dt.clear();
		if(!(x>>64)){u64 y=x;do{t.dt.push_back((y%E18));y/=E18;}while(y);}
		else do{t.dt.push_back(divbase(x));}while(x);
	}
	static Bint from128(i128 v){Bint t{};t.sig=v<0;put(t,v<0?0-static_cast<u128>(v):static_cast<u128>(v));return t;}
	static bool smallu(const Bint&x,u64&v){
		if(x.dt.size()!=1)return 0;
		v=x.dt[0];
		return 1;
	}
	static bool smallv(const Bint&x,i128&v){
		u64 u{};
		if(!smallu(x,u))return 0;
		v=x.sig?-static_cast<i128>(u):static_cast<i128>(u);
		return 1;
	}
	static bool abs_less(const Bint&x,const Bint&y){
		if(x.dt.size()!=y.dt.size())return x.dt.size()<y.dt.size();
		for(Z i=x.dt.size();i--;)if(x.dt[i]!=y.dt[i])return x.dt[i]<y.dt[i];
		return false;
	}
	static Bint abs_add(const Bint&a,const Bint&b){
		const V<T>&x{a.dt.size()>=b.dt.size()?a.dt:b.dt};
		const V<T>&y{a.dt.size()>=b.dt.size()?b.dt:a.dt};
		Bint t{};
		Z n{x.size()},m{y.size()};
		t.dt.resize(n+1);
		u64 c{};
		for(Z i{0};i<m;i++){u64 s{c+static_cast<u64>(x[i])+static_cast<u64>(y[i])};c=s/E18;t.dt[i]=(s-c*E18);}
		for(Z i{m};i<n;i++){u64 s{c+static_cast<u64>(x[i])};c=s/E18;t.dt[i]=(s-c*E18);}
		t.dt[n]=c;
		t.reload();
		return t;
	}
	static Bint abs_sub(const Bint&a,const Bint&b){
		Bint t{};
		Z n{a.dt.size()},m{b.dt.size()};
		t.dt.resize(n);
		u64 c{};
		Z i{};
		for(;i<m&&i<n;i++){long long s{static_cast<long long>(a.dt[i])-static_cast<long long>(b.dt[i])-static_cast<long long>(c)};if(s<0)s+=E18,c=1;else c=0;t.dt[i]=s;}
		for(;i<n;i++){long long s{static_cast<long long>(a.dt[i])-static_cast<long long>(c)};if(s<0)s+=E18,c=1;else c=0;t.dt[i]=s;}
		t.reload();
		return t;
	}
	static Bint mul_small(const Bint&a,i128 c){Bint t{};t.dt.resize(a.dt.size()+1);i128 k{};for(Z i{0};i<a.dt.size();i++){i128 v{a.dt[i]*c+k};t.dt[i]=(v%base);k=v/base;}t.dt[a.dt.size()]=k;t.reload();return t;}
	static void cv(const T*a,Z n,const T*b,Z m,u128*c){
		for(Z i{0};i<n;i++){u128 x=a[i];for(Z j{0};j<m;j++)c[i+j]+=x*b[j];}
	}
	static void nz(const u128*c,Z n,T*r){
		u128 carry{};
		for(Z i{0};i<n;i++){u128 s{c[i]+carry};r[i]=divbase(s);carry=s;}
	}
	static V<T> sm(const T*a,Z n,const T*b,Z m){
		Z L{n>m?n:m};
		V<T>r(L+1);
		u64 c{};
		for(Z i{0};i<L;i++){
			u64 s{c};
			if(i<n)s+=static_cast<u64>(a[i]);
			if(i<m)s+=static_cast<u64>(b[i]);
			c=s/E18;
			r[i]=(s-c*E18);
		}
		r[L]=c;
		while(r.size()>1&&!r.back())r.pop_back();
		return r;
	}
	static void sb(V<T>&x,const V<T>&y){
		i128 c{};
		for(Z i{0};i<x.size();i++){
			i128 s{static_cast<i128>(x[i])-(i<y.size()?static_cast<i128>(y[i]):0)-c};
			if(s<0)s+=E18,c=1;else c=0;
			x[i]=s;
		}
		while(x.size()>1&&!x.back())x.pop_back();
	}
	static void as(V<T>&r,const V<T>&z,Z k){
		if(r.size()<z.size()+k+1)r.resize(z.size()+k+1);
		u64 c{};
		Z i{0};
		for(;i<z.size();i++){u64 s{static_cast<u64>(r[k+i])+static_cast<u64>(z[i])+c};c=s/E18;r[k+i]=(s-c*E18);}
		for(Z j{k+i};c;j++){u64 s{static_cast<u64>(r[j])+c};c=s/E18;r[j]=(s-c*E18);}
	}
	static V<T> mv(const T*a,Z n,const T*b,Z m){
		if(n<m){const T*p{a};a=b;b=p;Z q{n};n=m;m=q;}
		V<T>r{};
		if(m<=320){
			V<u128>c(n+m,0);
			cv(a,n,b,m,c.data());
			r.resize(n+m);
			nz(c.data(),n+m,r.data());
			while(r.size()>1&&!r.back())r.pop_back();
			return r;
		}
		Z k{m>>1};
		V<T>z0{mv(a,k,b,k)},z2{mv(a+k,n-k,b+k,m-k)};
		V<T>s{sm(a,k,a+k,n-k)},t{sm(b,k,b+k,m-k)};
		V<T>z1{mv(s.data(),s.size(),t.data(),t.size())};
		sb(z1,z0);
		sb(z1,z2);
		r=z0;
		as(r,z1,k);
		as(r,z2,2*k);
		while(r.size()>1&&!r.back())r.pop_back();
		return r;
	}
	static constexpr u64 q1{998244353},q2{1004535809},q3{469762049},q12{1002772198720536577};
	static constexpr u64 n1{998244351},n2{1004535807},n3{469762047},u1{932051910},u2{542374313},u3{460175152},w12{669690699},w3{354521948};
	static constexpr Z th{2800};
	static u64 mm(u64 a,u64 b,u64 p,u64 np){
		u64 t{a*b},m{(t&4294967295ull)*np&4294967295ull},r{(t+m*p)>>32};
		return r>=p?r-p:r;
	}
	static u64 mp(u64 a,u64 e,u64 p,u64 np,u64 r2){
		u64 t{mm(a,r2,p,np)},r{mm(1,r2,p,np)};
		for(;e;e>>=1){if(e&1)r=mm(r,t,p,np);t=mm(t,t,p,np);}
		return r;
	}
	static void rt(V<u64>&t,Z n,u64 p,u64 np,u64 r2,bool iv){
		t.resize(n);
		for(Z l{2};l<=n;l<<=1){
			u64 e{iv?p-1-(p-1)/l:(p-1)/l},wl{mp(3,e,p,np,r2)},w{mm(1,r2,p,np)};
			for(Z j{0};j<l/2;j++){t[l/2+j]=w;w=mm(w,wl,p,np);}
		}
	}
	static void tf(V<u64>&f,Z n,const V<u64>&t,u64 p,u64 np,u64 r2,bool iv){
		for(Z i{1},j{};i<n;i++){
			Z b{n>>1};
			for(;j&b;b>>=1)j^=b;
			j^=b;
			if(i<j){u64 q{f[i]};f[i]=f[j];f[j]=q;}
		}
		for(Z l{2};l<=n;l<<=1)for(Z i{0};i<n;i+=l)for(Z j{0};j<l/2;j++){
			u64 u{f[i+j]},v{mm(f[i+j+l/2],t[l/2+j],p,np)},s{u+v};
			f[i+j]=s>=p?s-p:s;
			f[i+j+l/2]=u>=v?u-v:u+p-v;
		}
		if(iv){
			u64 s{mp(n%p,p-2,p,np,r2)};
			for(Z i{0};i<n;i++)f[i]=mm(f[i],s,p,np);
		}
	}
	static void nv(V<u64>&f,V<u64>&g,const u64*a,Z n,const u64*b,Z m,Z N,u64 p,u64 np,u64 r2,V<u64>&o){
		f.resize(N);g.resize(N);
		for(Z i{0};i<n;i++){u64 v{a[i]};while(v>=p)v-=p;f[i]=mm(v,r2,p,np);}
		for(Z i{n};i<N;i++)f[i]=0;
		for(Z i{0};i<m;i++){u64 v{b[i]};while(v>=p)v-=p;g[i]=mm(v,r2,p,np);}
		for(Z i{m};i<N;i++)g[i]=0;
		V<u64>t{};
		rt(t,N,p,np,r2,0);
		tf(f,N,t,p,np,r2,0);tf(g,N,t,p,np,r2,0);
		for(Z i{0};i<N;i++)f[i]=mm(f[i],g[i],p,np);
		rt(t,N,p,np,r2,1);
		tf(f,N,t,p,np,r2,1);
		o.resize(n+m-1);
		for(Z i{0};i<o.size();i++)o[i]=mm(f[i],1,p,np);
	}
	static V<T> nm(const Bint&a,const Bint&b){
		V<u64>ha(a.dt.size()*2),hb(b.dt.size()*2),f,g,r1,r2,r3;
		Z n{half(a,ha.data())},m{half(b,hb.data())},L{n+m-1},N{1};
		while(N<L)N<<=1;
		nv(f,g,ha.data(),n,hb.data(),m,N,q1,n1,u1,r1);
		nv(f,g,ha.data(),n,hb.data(),m,N,q2,n2,u2,r2);
		V<u64>c(L);
		for(Z i{0};i<L;i++)c[i]=static_cast<u64>(static_cast<u128>(r1[i])+q1*(static_cast<u128>(r2[i]>=r1[i]?r2[i]-r1[i]:r2[i]+q2-r1[i])*w12%q2));
		nv(f,g,ha.data(),n,hb.data(),m,N,q3,n3,u3,r3);
		V<u64>d{};
		d.reserve(L+4);
		u128 cy{};
		for(Z i{0};i<L;i++){
			u64 v{static_cast<u64>(static_cast<u128>(r3[i]>=c[i]%q3?r3[i]-c[i]%q3:r3[i]+q3-c[i]%q3)*w3%q3)};
			u128 s{static_cast<u128>(c[i])+static_cast<u128>(q12)*v+cy};
			d.push_back(div9(s));
			cy=s;
		}
		while(cy)d.push_back(div9(cy));
		V<T>r((d.size()+1)/2);
		for(Z i{0};i<r.size();i++)r[i]=(d[2*i]+(2*i+1<d.size()?d[2*i+1]*1000000000ull:0));
		while(r.size()>1&&!r.back())r.pop_back();
		return r;
	}
	friend Bint operator*(const Bint&a,const Bint&b){
		i128 x{},y{};
		if(smallv(a,x)&&smallv(b,y))return from128(x*y);
		Z n{a.dt.size()},m{b.dt.size()};
		Bint t{};
		t.sig=a.sig!=b.sig;
		if(n<=32&&m<=32){
			u128 c[64];
			for(Z i{0};i<n+m;i++)c[i]=0;
			cv(a.dt.data(),n,b.dt.data(),m,c);
			t.dt.resize(n+m);
			nz(c,n+m,t.dt.data());
		}else if((n<m?n:m)>=th&&n+m<=1048576)t.dt=nm(a,b);
		else t.dt=mv(a.dt.data(),n,b.dt.data(),m);
		t.reload();
		return t;
	}
	static Z half(const Bint&x,u64*h){
		Z k{};
		for(Z i{0};i<x.dt.size();i++){
			u64 v=x.dt[i];
			h[k++]=v%1000000000ull;
			h[k++]=v/1000000000ull;
		}
		while(k>1&&!h[k-1])k--;
		return k;
	}
	static constexpr Z bth=32;
	static u64 nf(u64 t){return t<500000000000000000ull?(500000000000000000ull+(t-1))/t:1;}
	static void tb(Bint&x){while(x.dt.size()>1&&!x.dt.back())x.dt.pop_back();if(x.dt.empty())x.dt.push_back({});}
	static Bint lo(const Bint&x,Z k){
		Bint t{};
		Z n{x.dt.size()<k?x.dt.size():k};
		if(n)t.dt.assign(x.dt.begin(),x.dt.begin()+n);
		return t;
	}
	static Bint hi(const Bint&x,Z k){
		Bint t{};
		if(x.dt.size()>k)t.dt.assign(x.dt.begin()+k,x.dt.end());
		return t;
	}
	static Bint cc(const Bint&h,const Bint&l,Z k){
		Bint t{};
		t.dt.assign(k,{});
		Z n{l.dt.size()<k?l.dt.size():k};
		for(Z i{0};i<n;i++)t.dt[i]=l.dt[i];
		for(Z i{0};i<h.dt.size();i++)t.dt.push_back(h.dt[i]);
		tb(t);
		return t;
	}
	static Bint bl(const Bint&x,Z o,Z k){
		Bint t{};
		Z b{o<x.dt.size()?o:x.dt.size()},e{o+k<x.dt.size()?o+k:x.dt.size()};
		if(b<e)t.dt.assign(x.dt.begin()+b,x.dt.begin()+e);
		return t;
	}
	static Bint ds(const Bint&x,u64 d){
		Bint t{};
		t.dt.resize(x.dt.size());
		u64 c{};
		for(Z i{x.dt.size()};i--;){
			u128 v{static_cast<u128>(c)*E18+x.dt[i]};
			u64 q=v/d;
			t.dt[i]=q;
			c=v-static_cast<u128>(q)*d;
		}
		tb(t);
		return t;
	}
	static void kd(const Bint&A,const Bint&C,Bint&q,Bint&r){
		Z n{A.dt.size()},m{C.dt.size()};
		if(m==1){
			u64 d=C.dt[0],k{};
			q.dt.resize(n);
			for(Z i{n};i--;){
				u128 x{static_cast<u128>(k)*E18+A.dt[i]};
				u64 t=x/d;
				q.dt[i]=t;
				k=x-static_cast<u128>(t)*d;
			}
			tb(q);
			r.dt.assign(1,k);
			return;
		}
		u64 d{nf(C.dt[m-1])},k{};
		V<u64>u(n+1),v(m);
		if(d==1){
			for(Z i{0};i<n;i++)u[i]=A.dt[i];
			for(Z i{0};i<m;i++)v[i]=C.dt[i];
		}else{
			for(Z i{0};i<n;i++){u128 p{static_cast<u128>(A.dt[i])*d+k};u[i]=p%E18;k=p/E18;}
			u[n]=k;k=0;
			for(Z i{0};i<m;i++){u128 p{static_cast<u128>(C.dt[i])*d+k};v[i]=p%E18;k=p/E18;}
		}
		q.dt.assign(n-m+1,0);
		for(Z j{n-m+1};j--;){
			u128 nm{static_cast<u128>(u[j+m])*E18+u[j+m-1]};
			u64 qh=nm/static_cast<u128>(v[m-1]),rh=nm%static_cast<u128>(v[m-1]);
			if(qh>=E18){qh=999999999999999999ull;rh=nm-static_cast<u128>(qh)*v[m-1];}
			while(rh<E18&&static_cast<u128>(qh)*v[m-2]>static_cast<u128>(rh)*E18+u[j+m-2])qh--,rh+=v[m-1];
			u64 c2{},br{};
			for(Z i{0};i<m;i++){
				u128 p{static_cast<u128>(qh)*v[i]+c2};
				c2=p/E18;
				u64 pl=p-c2*E18,t{u[j+i]-pl-br};
				br=t>>63;
				u[j+i]=t+br*E18;
			}
			u64 t0{u[j+m]-c2-br};
			br=t0>>63;
			u[j+m]=t0+br*E18;
			while(br){
				qh--;
				u64 g{};
				for(Z i{0};i<m;i++){u64 s{u[j+i]+v[i]+g};if(s>=E18){s-=E18;g=1;}else g=0;u[j+i]=s;}
				u[j+m]+=g;
				if(u[j+m]>=E18)u[j+m]-=E18;
				br=u[j+m]!=0;
			}
			q.dt[j]=qh;
		}
		tb(q);
		r.dt.resize(m);
		if(d==1){for(Z i{0};i<m;i++)r.dt[i]=u[i];}
		else{u64 g{};for(Z i{m};i--;){u128 x{static_cast<u128>(g)*E18+u[i]};u64 t=x/d;r.dt[i]=t;g=x-static_cast<u128>(t)*d;}}
		tb(r);
	}
	static void bz(const Bint&A,const Bint&C,Bint&q,Bint&r){
		Z n{A.dt.size()},m{C.dt.size()};
		if(n<m){q=Bint{};r=A;return;}
		if(m<=bth){kd(A,C,q,r);return;}
		Z h{m>>1},K{(n-m+1+h-1)/h};
		Bint R{bl(A,h*(K-1),n-h*(K-1))};
		q.dt.assign(h*K,0);
		for(Z i{K};i--;){
			if(i!=K-1)R=cc(R,bl(A,h*i,h),h);
			Bint qa{},ra{};
			if(R.dt.size()<=C.dt.size()){
				if(abs_less(R,C))qa=Bint{},ra=R;
				else qa=Bint{1},ra=R-C;
			}else{
				Bint b1{hi(C,h)},b0{lo(C,h)},ap{hi(R,h)},rl{lo(R,h)};
				bz(ap,b1,qa,ra);
				Bint t{cc(ra,rl,h)};
				t-=qa*b0;
				while(t.sig)qa-=Bint{1},t+=C;
				while(!abs_less(t,C))qa+=Bint{1},t-=C;
				ra=t;
			}
			for(Z j{0};j<qa.dt.size()&&h*i+j<q.dt.size();j++)q.dt[h*i+j]=qa.dt[j];
			R=ra;
		}
		tb(q);
		r=R;
	}
	static void divmod(const Bint&a,const Bint&b,Bint&q,Bint&r){
		Bint A{a},C{b};
		A.sig=0;C.sig=0;
		A.reload();C.reload();
		if(!A||abs_less(A,C)){
			q=Bint{};
			r=A;
			if(r)r.sig=a.sig;
			return;
		}
		Z m{C.dt.size()};
		if(m<=bth)kd(A,C,q,r);
		else{
			u64 d{nf(C.dt[m-1])};
			if(d!=1)A=mul_small(A,d),C=mul_small(C,d);
			bz(A,C,q,r);
			if(d!=1)r=ds(r,d);
		}
		q.reload();r.reload();
		if(q)q.sig=a.sig!=b.sig;
		if(r)r.sig=a.sig;
	}
	friend Bint operator+(const Bint&a,const Bint&b){
		i128 x{},y{};
		if(smallv(a,x)&&smallv(b,y))return from128(x+y);
		Bint t{};
		if(a.sig==b.sig)t=abs_add(a,b),t.sig=a.sig;
		else if(abs_less(a,b))t=abs_sub(b,a),t.sig=b.sig;
		else t=abs_sub(a,b),t.sig=a.sig;
		t.reload();
		return t;
	}
	Bint&operator+=(const Bint&x){*this=*this+x;return*this;}
	friend Bint operator-(const Bint&a){Bint t{a};if(t)t.sig=!t.sig;return t;}
	friend Bint operator-(const Bint&a,const Bint&b){
		i128 x{},y{};
		if(smallv(a,x)&&smallv(b,y))return from128(x-y);
		return a+(-b);
	}
	Bint&operator-=(const Bint&x){*this=*this-x;return*this;}
	Bint&operator*=(const Bint&x){*this=*this*x;return*this;}
	friend Bint operator/(const Bint&a,const Bint&b){Bint q{},r{};divmod(a,b,q,r);return q;}
	Bint&operator/=(const Bint&x){*this=*this/x;return*this;}
	friend Bint operator%(const Bint&a,const Bint&b){Bint q{},r{};divmod(a,b,q,r);return r;}
	Bint&operator%=(const Bint&x){*this=*this%x;return*this;}
	friend bool operator<(const Bint&a,const Bint&b){
		return a.sig!=b.sig?a.sig:(a.sig?abs_less(b,a):abs_less(a,b));
	}
	friend bool operator<=(const Bint&a,const Bint&b){return!(b<a);}
	friend bool operator==(const Bint&a,const Bint&b){return a.sig==b.sig&&a.dt==b.dt;}
	friend bool operator!=(const Bint&a,const Bint&b){return!(a==b);}
	friend bool operator>=(const Bint&a,const Bint&b){return!(a<b);}
	friend bool operator>(const Bint&a,const Bint&b){return b<a;}
};
using BI=Bint<i128>;
template<typename U>struct fracb_isbint{static constexpr bool value=false;};
template<typename U>struct fracb_isbint<Bint<U>>{static constexpr bool value=true;};
template<typename T=BI>
class Frac{
	public:
	T a,b;
	Frac():a{},b{1}{}
	Frac(const T&_a,const T&_b=1):a{_a},b{_b}{simplify();}
	template<typename U=T,typename=std::enable_if_t<fracb_isbint<U>::value>>
	Frac(const i128&x):a{x},b{1}{}
	static constexpr __float128 eps=fpe<T>::v;
	static bool one(const T&x){
		if constexpr(fracb_isbint<T>::value)return !x.sig&&x.dt.size()==1&&x.dt[0]==1;
		return x==1;
	}
	static u128 wd(const T&x){return x.dt.size()>1?static_cast<u128>(x.dt[1])*static_cast<u128>(T::base)+static_cast<u128>(x.dt[0]):static_cast<u128>(x.dt[0]);}
	static T mkw(u128 v){
		T t{};
		t.dt.clear();
		do{t.dt.push_back(v%static_cast<u128>(T::base));v/=static_cast<u128>(T::base);}while(v);
		return t;
	}
	static u128 gw(u128 x,u128 y){
		while(y){u128 t{x%y};x=y;y=t;}
		return x;
	}
	static T ed(const T&a,u128 m){
		T t{a};u128 c{0};
		for(Z i=t.dt.size();i--;){u128 v{c*static_cast<u128>(T::base)+static_cast<u128>(t.dt[i])};t.dt[i]=(v/m);c=v%m;}
		while(t.dt.size()>1&&!t.dt.back())t.dt.pop_back();
		return t;
	}
	static T dv(const T&a,const T&g){
		if(a.dt==g.dt)return a.sig?T{-1}:T{1};
		if(g.dt.size()==1)return ed(a,(u128)g.dt[0]);
		if(a.dt.size()<=2&&g.dt.size()<=2){T q{mkw(wd(a)/wd(g))};q.sig=a.sig;return q;}
		return a/g;
	}
	static void hl(T&x){i128 c{0};
		for(Z i=x.dt.size();i--;){i128 v{static_cast<i128>(x.dt[i])+c*T::base};x.dt[i]=v/2;c=v&1;}
		while(x.dt.size()>1&&!x.dt.back())x.dt.pop_back();
	}
	static void dl(T&x){
		i128 c{0};
		for(Z i=0;i<x.dt.size();i++){i128 v{static_cast<i128>(x.dt[i])*2+c};x.dt[i]=v%T::base;c=v/T::base;}
		if(c)x.dt.push_back(c);
	}
	static void sb(T&x,const T&y){i128 c{0};
		for(Z i=0;i<x.dt.size();i++){i128 v{static_cast<i128>(x.dt[i])-c-(i<y.dt.size()?static_cast<i128>(y.dt[i]):0)};if(v<0)v+=T::base,c=1;else c=0;x.dt[i]=v;}
		while(x.dt.size()>1&&!x.dt.back())x.dt.pop_back();
	}
	static T bg2(T x,T y){
		Z z{0},w{0};
		while(!(x.dt[0]&1)){hl(x);z++;}
		while(!(y.dt[0]&1)){hl(y);w++;}
		if(z>w)z=w;
		for(;;){if(T::abs_less(x,y))std::swap(x,y);if(x.dt==y.dt)break;sb(x,y);while(!(x.dt[0]&1))hl(x);}
		while(z--)dl(x);return x;
	}
	static T gcd(T x,T y){
		if constexpr(fpt<T>){
			if(x<0)x=-x;
			if(y<0)y=-y;
			while(y){T t{fpf::mod(x,y)};x=y;y=t;}
			return x;
		}else if constexpr(fracb_isbint<T>::value){
			for(;;){
				if(T::abs_less(x,y))std::swap(x,y);
				if(!y)return x;
				if(y.dt.size()==1){
					u128 m{static_cast<u128>(y.dt[0])},r{0};if(m<2)return mkw(1);
					for(Z i=x.dt.size();i--;)r=(r*static_cast<u128>(T::base)+static_cast<u128>(x.dt[i]))%m;
					return mkw(gw(r,m));
				}
				if(x.dt.size()<=2)return mkw(gw(wd(x),wd(y)));
				if(x.dt.size()<=4)return bg2(x,y);
				T r{x%y};x=y;y=r;
			}
		}else{
			if(x<0)x=-x;
			if(y<0)y=-y;
			while(y){T t{x%y};x=y;y=t;}
			return x;
		}
	}
	static bool eq(const T&x,const T&y){
		if constexpr(is_same_v<T,__float128>)return fpf::eq(x,y,eps);
		if constexpr(fpt<T>)return fpf::abs(x-y)<eps;
		return x==y;
	}
	void simplify(){
		if(b<0)a=-a,b=-b;
		if constexpr(fpt<T>)return;
		else{
			if(one(b))return;
			T g{gcd(a,b)};
			if(g>1)a=dv(a,g),b=dv(b,g);
		}
	}
	explicit operator bool()const{return !!a;}
	static Frac qt(const T&p,const T&q,const T&r,const T&s){
		T g{gcd(p,s)};
		if(g>1){T h{gcd(r,q)};if(h>1)return{dv(p,g)*dv(r,h),dv(q,h)*dv(s,g)};return{dv(p,g)*r,q*dv(s,g)};}
		T h{gcd(r,q)};
		if(h>1)return{p*dv(r,h),dv(q,h)*s};
		return{p*r,q*s};
	}
	static Frac ad(const Frac&x,const Frac&y,bool m){
		if constexpr(fracb_isbint<T>::value){
			if(x.b.dt.size()+y.b.dt.size()>=3&&x.b&&y.b){
				T g{gcd(x.b,y.b)};
				if(g>1){T g1{dv(x.b,g)},g2{dv(y.b,g)};return{m?x.a*g2-y.a*g1:x.a*g2+y.a*g1,g1*y.b};}
			}
		}
		return{m?x.a*y.b-y.a*x.b:x.a*y.b+y.a*x.b,x.b*y.b};
	}
	friend Frac operator+(const Frac&x,const Frac&y){return ad(x,y,0);}
	friend Frac operator-(const Frac&x,const Frac&y){return ad(x,y,1);}
	friend Frac operator-(const Frac&x){return{-x.a,x.b};}
	friend Frac operator*(const Frac&x,const Frac&y){
		if constexpr(fracb_isbint<T>::value)if(x.a.dt.size()+y.a.dt.size()>=5)return qt(x.a,x.b,y.a,y.b);
		return{x.a*y.a,x.b*y.b};
	}
	friend Frac operator/(const Frac&x,const Frac&y){
		if constexpr(fracb_isbint<T>::value)if(x.a.dt.size()+y.b.dt.size()>=5)return qt(x.a,x.b,y.b,y.a);
		return{x.a*y.b,x.b*y.a};
	}
	Frac&operator+=(const Frac&x){*this=*this+x;return *this;}
	Frac&operator-=(const Frac&x){*this=*this-x;return *this;}
	Frac&operator*=(const Frac&x){*this=*this*x;return *this;}
	Frac&operator/=(const Frac&x){*this=*this/x;return *this;}
	friend bool operator<(const Frac&x,const Frac&y){return x.a*y.b<y.a*x.b;}
	friend bool operator<=(const Frac&x,const Frac&y){return !(y<x);}
	friend bool operator==(const Frac&x,const Frac&y){return eq(x.a*y.b,y.a*x.b);}
	friend bool operator!=(const Frac&x,const Frac&y){return !(x==y);}
	friend bool operator>=(const Frac&x,const Frac&y){return !(x<y);}
	friend bool operator>(const Frac&x,const Frac&y){return y<x;}
};
FPTP struct pex{using type=T;static constexpr bool v=false;};
template<typename U>struct pex<Frac<Bint<U>>>{using type=Bint<U>;static constexpr bool v=true;};
template<typename T=Frac<BI>>
class Poly{
	public:
	V<T>dt;
	void assign(Z _n,const T&x){dt.assign(_n,x);}
	T&operator[](Z i){return dt[i];}
	const T&operator[](Z i)const{return dt[i];}
	auto begin(){return dt.begin();}
	auto end(){return dt.end();}
	auto begin()const{return dt.begin();}
	auto end()const{return dt.end();}
	Poly(){}
	Poly(const V<T>&_dt):dt{_dt}{}
	Poly(std::initializer_list<T>_dt):dt{_dt}{}
	Poly(const T&t):dt{t}{}
	static constexpr __float128 eps=fpe<T>::v;
	static bool zero(const T&x){if constexpr(cxt<T>::value)return fpf::cabs(x)<eps;else if constexpr(fpt<T>)return fpf::abs(x)<eps;else return !x;}
	static bool eq(const T&a,const T&b){if constexpr(cxt<T>::value)return fpf::cabs(a-b)<eps;if constexpr(is_same_v<T,__float128>)return fpf::eq(a,b,eps);if constexpr(fpt<T>)return fpf::abs(a-b)<eps;return a==b;}
	void reload(){while(!dt.empty()&&zero(dt.back()))dt.pop_back();if(dt.empty())dt.push_back({});}
	T operator()(T x)const{T r{};for(Z i=dt.size();i--;)r=r*x+dt[i];return r;}
	template<typename A>
	static A ig(A x,A y){if(x.sig)x=-x;if(y.sig)y=-y;while(y){A t{x%y};x=y;y=t;}return x;}
	template<typename A>
	static A lc(const V<T>&v){A l{1};for(const T&x:v){if(x.b==l||(x.b.dt.size()==1&&x.b.dt[0]==1))continue;l=l/ig(l,x.b)*x.b;if(l.dt.size()>8)return A{};}return l;}
	template<typename A>
	static V<A> conv(const V<A>&a,const V<A>&b){
		V<A>u{a},v{b},c;
		if(u.size()<v.size())u.swap(v);
		c.resize(u.size()+v.size()-1);
		if(v.size()<=32){for(Z i=0;i<u.size();i++)for(Z j=0;j<v.size();j++)c[i+j]+=u[i]*v[j];return c;}
		Z h{v.size()/2};
		V<A>u1(u.begin()+h,u.end()),v1(v.begin()+h,v.end()),z0,z1,z2;
		u.resize(h);v.resize(h);
		z0=conv(u,v);z2=conv(u1,v1);
		u.resize(u1.size());v.resize(v1.size());
		for(Z i=0;i<u1.size();i++)u[i]+=u1[i];for(Z i=0;i<v1.size();i++)v[i]+=v1[i];
		z1=conv(u,v);
		for(Z i=0;i<z0.size();i++)z1[i]-=z0[i];for(Z i=0;i<z2.size();i++)z1[i]-=z2[i];
		for(Z i=0;i<z1.size();i++){if(i<z0.size())c[i]+=z0[i];c[i+h]+=z1[i];c[i+2*h]+=z2[i];}
		return c;
	}
	static bool smul(const Poly&a,const Poly&b,Poly&t){
		using A=typename pex<T>::type;
		A la{lc<A>(a.dt)},lb{lc<A>(b.dt)};
		if(!la||!lb)return false;
		V<A>ua(a.dt.size()),ub(b.dt.size());
		for(Z i=0;i<ua.size();i++)ua[i]=a.dt[i].b==la?a.dt[i].a:a.dt[i].a*(la/a.dt[i].b);
		for(Z i=0;i<ub.size();i++)ub[i]=b.dt[i].b==lb?b.dt[i].a:b.dt[i].a*(lb/b.dt[i].b);
		V<A>uc{conv(ua,ub)};
		A den{la*lb};
		t.dt.resize(uc.size());
		if(den.dt.size()==1&&den.dt[0]==1&&!den.sig)for(Z i=0;i<uc.size();i++)t.dt[i].a=uc[i];
		else for(Z i=0;i<uc.size();i++)t.dt[i]=T{uc[i],den};
		t.reload();
		return true;
	}
	static void pdiv(const V<T>&a,const V<T>&b,V<T>&q){
		Z m{b.size()},sz{a.size()};
		V<T>r{a};
		q.assign(sz,T{});
		while(sz>=m&&(sz>1||!zero(r[0]))){
			Z sh{sz-m};
			T c{r[sz-1]/b[m-1]};
			q[sh]=c;
			if constexpr(pex<T>::v){
				for(Z i=0;i+1<m;i++){
					T&ri{r[sh+i]};
					auto cb{c.b*b[i].b};
					ri=T{ri.a*cb-c.a*b[i].a*ri.b,ri.b*cb};
				}
				r[sh+m-1]=T{};
			}else for(Z i=0;i<m;i++)r[sh+i]-=c*b[i];
			while(sz>1&&zero(r[sz-1]))sz--;
		}
	}
	friend Poly operator+(const Poly&a,const Poly&b){
		Poly t{};
		t.dt.resize(std::max(a.dt.size(),b.dt.size()));
		for(Z i=0;i<t.dt.size();i++){if(i<a.dt.size())t.dt[i]+=a.dt[i];if(i<b.dt.size())t.dt[i]+=b.dt[i];}
		t.reload();
		return t;
	}
	Poly&operator+=(const Poly&b){
		if(dt.size()<b.dt.size())dt.resize(b.dt.size());
		for(Z i=0;i<b.dt.size();i++)dt[i]+=b.dt[i];
		reload();
		return *this;
	}
	friend Poly operator-(const Poly&a){
		Poly t{a};
		for(T&x:t.dt)x=-x;
		t.reload();
		return t;
	}
	friend Poly operator-(const Poly&a,const Poly&b){Poly t{a};t-=b;return t;}
	Poly&operator-=(const Poly&b){
		if(dt.size()<b.dt.size())dt.resize(b.dt.size());
		for(Z i=0;i<b.dt.size();i++)dt[i]-=b.dt[i];
		reload();
		return *this;
	}
	friend Poly operator*(const Poly&a,const Poly&b){
		Poly t{};
		if constexpr(pex<T>::v)if(!a.dt.empty()&&!b.dt.empty()&&smul(a,b,t))return t;
		t.dt.resize(a.dt.size()+b.dt.size()-1);
		for(Z i=0;i<a.dt.size();i++)for(Z j=0;j<b.dt.size();j++)t.dt[i+j]+=a.dt[i]*b.dt[j];
		t.reload();
		return t;
	}
	Poly&operator*=(const Poly&b){*this=*this*b;return *this;}
	friend Poly operator/(const Poly&a,const Poly&b){Poly t{a};t/=b;return t;}
	Poly&operator/=(const Poly&b){Poly t{};pdiv(dt,b.dt,t.dt);dt.swap(t.dt);reload();return *this;}
	friend Poly operator%(const Poly&a,const Poly&b){Poly t{a};t%=b;return t;}
	Poly&operator%=(const Poly&b){*this-=*this/b*b;reload();return *this;}
	friend bool operator<(const Poly&a,const Poly&b){
		Poly x{a},y{b};
		x.reload();y.reload();
		if(x.dt.size()!=y.dt.size())return x.dt.size()<y.dt.size();
		for(Z i=x.dt.size();i--;)if(!eq(x.dt[i],y.dt[i]))return x.dt[i]<y.dt[i];
		return false;
	}
	friend bool operator<=(const Poly&a,const Poly&b){return !(b<a);}
	friend bool operator==(const Poly&a,const Poly&b){
		Poly x{a},y{b};
		x.reload();y.reload();
		if(x.dt.size()!=y.dt.size())return false;
		for(Z i=0;i<x.dt.size();i++)if(!eq(x.dt[i],y.dt[i]))return false;
		return true;
	}
	friend bool operator!=(const Poly&a,const Poly&b){return !(a==b);}
	friend bool operator>=(const Poly&a,const Poly&b){return !(a<b);}
	friend bool operator>(const Poly&a,const Poly&b){return b<a;}
};
namespace fpg{
using FB=BI;
using FF=Frac<FB>;
using FP=Poly<FF>;
static const u32 pf[6]={1000000007u,1000000009u,998244353u,1004535809u,469762049u,167772161u};
inline u64 pw(u64 a,u64 e,u64 p){u64 r{1};a%=p;while(e){if(e&1)r=r*a%p;a=a*a%p;e>>=1;}return r;}
inline u32 mb(const FB&x,u32 p){
	u64 r{},b=FB::base%p;
	for(Z i=x.dt.size();i--;)r=(r*b+(x.dt[i]%p))%p;
	if(x.sig&&r)r=p-r;
	return r;
}
inline bool red(const FP&a,u32 p,V<u32>&v){
	v.resize(a.dt.size());
	for(Z i=0;i<a.dt.size();i++){
		u32 n=mb(a.dt[i].a,p),d=mb(a.dt[i].b,p);
		if(!d)return 0;
		v[i]=(n*pw(d,p-2,p)%p);
	}
	while(v.size()>1&&!v.back())v.pop_back();
	return v.size()>1||v[0];
}
inline void gcdp(V<u32>a,V<u32>b,u32 p,V<u32>&g){
	while(!(b.size()==1&&!b[0])){
		V<u32>r{a};
		u32 v=pw(b.back(),p-2,p);
		while(r.size()>=b.size()&&!(r.size()==1&&!r[0])){
			Z sh=r.size()-b.size();
			u32 c=(static_cast<u64>(r.back())*v%p);
			if(c)for(Z i=0;i<b.size();i++){
				u32 t=(static_cast<u64>(c)*b[i]%p);
				r[sh+i]=r[sh+i]>=t?(r[sh+i]-t):(r[sh+i]+p-t);
			}
			while(r.size()>1&&!r.back())r.pop_back();
		}
		a=b;b=r;
	}
	u32 v=pw(a.back(),p-2,p);
	for(auto&x:a)x=(static_cast<u64>(x)*v%p);
	g=a;
}
inline FP monic(FP x){x.reload();if(x.dt.size()>1||!FP::zero(x.dt[0]))return x/x.dt.back();return x;}
inline bool fast(const FP&aa,const FP&bb,FP&o){
	FP a{aa},b{bb};
	a.reload();b.reload();
	if(b.dt.size()==1&&FP::zero(b.dt[0])){o=monic(a);return 1;}
	if(a.dt.size()==1&&FP::zero(a.dt[0])){o=monic(b);return 1;}
	if(a.dt.size()==1||b.dt.size()==1){o=FP{1};return 1;}
	FB la{a.dt.back().a},lb{b.dt.back().a},M{1};
	V<FB>X{};
	int k{-1};
	for(int it=0;it<6;it++){
		u32 p=pf[it];
		if(!mb(la,p)||!mb(lb,p))continue;
		V<u32>A,B,g;
		if(!red(a,p,A)||!red(b,p,B))continue;
		if(A.size()!=a.dt.size()||B.size()!=b.dt.size())continue;
		gcdp(A,B,p,g);
		int kk=g.size()-1;
		if(!kk){o=FP{1};return 1;}
		if(k<0)k=kk;
		if(kk!=k)continue;
		u32 c=mb(la,p);
		V<u32>r(g.size());
		for(Z j=0;j<g.size();j++)r[j]=(static_cast<u64>(c)*g[j]%p);
		if(X.empty()){
			X.resize(g.size());
			for(Z j=0;j<g.size();j++)X[j]=FB{r[j]};
			M=FB{p};
		}else{
			u32 iv=pw(mb(M,p),p-2,p);
			for(Z j=0;j<g.size();j++){
				u32 x=mb(X[j],p),t=x<=r[j]?(r[j]-x):(r[j]+p-x);
				u32 u=(static_cast<u64>(t)*iv%p);
				if(u)X[j]+=M*FB{u};
			}
			M*=FB{p};
		}
		if(it){
			FP q{};
			q.dt.resize(X.size());
			for(Z j=0;j<X.size();j++){
				FB v{X[j]};
				if(v*FB{2}>M)v-=M;
				q.dt[j]=FF{v,la};
			}
			q.reload();
			if(q.dt.size()==X.size()){
				FP ra{a%q},rb{b%q};
				if((ra.dt.size()==1&&FP::zero(ra.dt[0]))&&(rb.dt.size()==1&&FP::zero(rb.dt[0]))){o=monic(q);return 1;}
			}
		}
	}
	return 0;
}
}
namespace std{
	FPTP
	Poly<T>gcd(const Poly<T>&a,const Poly<T>&b){
		if constexpr(is_same_v<T,Frac<BI>>){Poly<T>r{};if(fpg::fast(a,b,r))return r;}
		Poly<T>x{a},y{b};
		x.reload();y.reload();
		while(y.dt.size()>1||!Poly<T>::zero(y.dt[0])){Poly<T>r{x%y};x=move(y);y=move(r);}
		if(x.dt.size()>1||!Poly<T>::zero(x.dt[0]))x/=Poly<T>{x.dt.back()};
		return x;
	}
}
FPTP
Poly<T>Deri(Poly<T>p){
	if(p.dt.size()==1){p.dt[0]={};return p;}
	if constexpr(pex<T>::v){
		using A=typename pex<T>::type;
		for(Z i=1;i<p.dt.size();i++)p.dt[i-1]=T{p.dt[i].a*A(i),p.dt[i].b};
	}else for(Z i=1;i<p.dt.size();i++)p.dt[i-1]=p.dt[i]*T(i);
	p.dt.pop_back();
	p.reload();
	return p;
}
FPTP
Poly<T>inter(Poly<T>p){
	p.dt.push_back({});
	if constexpr(pex<T>::v){
		using A=typename pex<T>::type;
		for(Z i=p.dt.size()-1;i--;)p.dt[i+1]=T{p.dt[i].a,p.dt[i].b*A(i+1)};
	}else for(Z i=p.dt.size()-1;i--;)p.dt[i+1]=p.dt[i]/T(i+1);
	p.dt[0]={};
	p.reload();
	return p;
}
FPTP
Poly<T>radic(const Poly<T>&a){
	if(a.dt.size()==1&&Poly<T>::zero(a.dt[0]))return a;
	Poly<T>r{a/std::gcd(a,Deri(a))};
	if(r.dt.size()==1&&Poly<T>::zero(r.dt[0]))return r;
	return r/r.dt.back();
}
template<typename T1,typename T2>
T1 fpow(T1 a,T2 b){
	T1 ans{1};
	while(b){
		if(b&1)ans*=a;
		b>>=1;
		if(b)a*=a;
	}
	return ans;
}
template<typename T1,typename T2>
T1 fpow(T1 a,T2 b,T1 p){
	T1 ans{1};
	if constexpr(fpt<T1>){
		a=fpf::mod(a,p);
		while(b){
			if(b&1)ans=fpf::mod(ans*a,p);
			b>>=1;
			if(b)a=fpf::mod(a*a,p);
		}
	}else{
		a%=p;
		while(b){
			if(b&1)(ans*=a)%=p;
			b>>=1;
			if(b)(a*=a)%=p;
		}
	}
	return ans;
}
namespace dfpio{
#ifdef _getchar_nolock
	inline int rdc(){return _getchar_nolock();}
#else
	inline int rdc(){return getchar();}
#endif
	inline void dig(V<i128>&v,Z&n,int d){
		if(n%18==0)v.push_back(d);
		else v.back()=v.back()*10+d;
		n++;
	}
	FPTP
	inline void bld(const V<i128>&v,Z n,Bint<T>&t,bool neg){
		if(!n){t.dt.assign(1,{});t.sig=0;return;}
		Z r{n%18};if(!r)r=18;
		i128 cur{v.back()},p{1};
		for(Z j=0;j<r;j++)p*=10;
		t.dt.clear();t.dt.reserve(v.size()+1);
		for(Z i=v.size()-1;i--;){
			cur+=v[i]*p;
			t.dt.push_back((cur%E18));
			cur/=E18;
		}
		t.dt.push_back(cur);
		t.sig=neg;
		t.reload();
	}
	inline char rdg(V<i128>&v,Z&n,char c){
		while(c!=EOF&&!(c<'0'||c>'9')){dig(v,n,c-'0');c=rdc();}
		return c;
	}
	inline char rds(string&s,char c){
		while(c!=EOF&&!(c<'0'||c>'9')){s+=c;c=rdc();}
		return c;
	}
	inline void wr(const char*b,Z n){fwrite(b,1,n,stdout);}
	inline void puti(char*b,Z&p,u64 m,int w){
		static const char t[]="00010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899";
		char d[20];int n{0};
		while(m>=100){u32 r=m%100;m/=100;d[n++]=t[r*2+1];d[n++]=t[r*2];}
		if(m>9){d[n++]=t[m*2+1];d[n++]=t[m*2];}else d[n++]='0'+m;
		int q{w-n};if(q<0)q=0;
		while(q--)b[p++]='0';
		while(n)b[p++]=d[--n];
	}
	inline int bdig(const BI&x){
		u64 v=x.dt.back();
		int n=(x.dt.size()-1)*18;
		do{n++;v/=10;}while(v);
		return n;
	}
	inline BI bpow(BI a,int e){
		BI r{1};
		while(e){
			if(e&1)r*=a;
			e>>=1;
			if(e)a*=a;
		}
		return r;
	}
	inline void qwr(const __float128&x){
		u64 u[2];
		__builtin_memcpy(u,&x,16);
		bool neg{(u[1]>>63)!=0};
		int ex=((u[1]>>48)&32767);
		u128 m{(static_cast<u128>(u[1]&0xffffffffffffull)<<64)|u[0]};
		if(ex==32767){if(neg)putchar('-');fputs(m?"nan":"inf",stdout);return;}
		if(!ex&&!m){if(neg)putchar('-');putchar('0');return;}
		if(ex)m|=(u128)1<<112;
		int b{ex?ex-16495:-16494};
		BI N{static_cast<i128>(m)},D{1};
		if(b>=0)N*=bpow(BI{2},b);
		else D=bpow(BI{2},-b);
		int nd{bdig(N)},dd{bdig(D)},E0{nd-dd},E;
		bool ge;
		if(E0>=0){
			BI P{bpow(BI{10},E0)};
			ge=!(N<D*P);
		}else{
			BI P{bpow(BI{10},-E0)};
			ge=!(N*P<D);
		}
		E=ge?E0:E0-1;
		int s{35-E};
		BI Ns{N},Ds{D};
		if(s>=0)Ns*=bpow(BI{10},s);
		else Ds*=bpow(BI{10},-s);
		BI Q{},R{};
		BI::divmod(Ns,Ds,Q,R);
		BI R2{BI::mul_small(R,2)};
		if(!(R2<Ds))Q+=BI{1};
		char c[64];Z p{};
		puti(c,p,Q.dt.back(),0);
		for(Z i=Q.dt.size()-1;i--;)puti(c,p,Q.dt[i],18);
		if(p>36){p=36;E++;}
		char d[40];
		int lead=36-p;
		for(int i=0;i<36;i++)d[i]=i<lead?'0':c[i-lead];
		int L{36};while(L>1&&d[L-1]=='0')L--;
		if(neg)putchar('-');
		if(E>=0&&E<36){
			if(E<L){for(int i=0;i<=E;i++)putchar(d[i]);if(E+1<L){putchar('.');for(int i=E+1;i<L;i++)putchar(d[i]);}}
			else{for(int i=0;i<L;i++)putchar(d[i]);for(int i=L;i<=E;i++)putchar('0');}
		}else if(E<0&&E>=-6){
			int z{-E-1};fputs("0.",stdout);for(int i=0;i<z;i++)putchar('0');for(int i=0;i<L;i++)putchar(d[i]);
		}else{
			putchar(d[0]);if(L>1){putchar('.');for(int i=1;i<L;i++)putchar(d[i]);}putchar('e');printf("%d",E);
		}
	}
	inline __float128 qmk(bool neg,int ex,u128 fr){
		u64 u[2]{static_cast<u64>(fr),(static_cast<u64>(neg)<<63)|(static_cast<u64>(ex)<<48)|static_cast<u64>(fr>>64)};
		__float128 r;__builtin_memcpy(&r,u,16);return r;
	}
	inline __float128 qdec(bool neg,u128 u,int q){
		if(!u)return qmk(neg,0,0);
		int du{1};u128 w{u};while(w>=10){w/=10;du++;}
		int E0{du+q-1};
		if(E0>4933)return qmk(neg,32767,0);
		if(E0<-4967)return qmk(neg,0,0);
		BI U{static_cast<i128>(u)},D{1};
		if(q>=0)U*=bpow(BI{10},q);
		else D=bpow(BI{10},-q);
		int e2=std::floor(static_cast<long double>(E0)*3.321928094887362347870319429489390L);
		for(;;){
			if(e2>=0){
				BI P{bpow(BI{2},e2)},T{D*P};
				if(U<T){e2--;continue;}
				BI T2{BI::mul_small(T,2)};
				if(!(U<T2)){e2++;continue;}
				break;
			}else{
				BI P{bpow(BI{2},-e2)},T{U*P};
				if(T<D){e2--;continue;}
				BI T2{BI::mul_small(D,2)};
				if(!(T<T2)){e2++;continue;}
				break;
			}
		}
		if(e2>16383)return qmk(neg,32767,0);
		if(e2>=-16382){
			int sh{112-e2};
			BI Ns{U},Ds{D};
			if(sh>=0)Ns*=bpow(BI{2},sh);
			else Ds*=bpow(BI{2},-sh);
			BI Q{},R{};
			BI::divmod(Ns,Ds,Q,R);
			u128 M{(u128)Q.dt[0]};
			if(Q.dt.size()>1)M+=(u128)Q.dt[1]*(u128)BI::base;
			BI R2{BI::mul_small(R,2)};
			int c=R2<Ds?0:(Ds<R2?1:(M&1));
			if(c)M++;
			if(M>>113){M>>=1;e2++;}
			int ex{e2+16383};
			if(ex>=32767)return qmk(neg,32767,0);
			return qmk(neg,ex,M&(((u128)1<<112)-1));
		}else{
			BI Ns{U*bpow(BI{2},16494)},Ds{D};
			BI Q{},R{};
			BI::divmod(Ns,Ds,Q,R);
			u128 M{(u128)Q.dt[0]};
			if(Q.dt.size()>1)M+=(u128)Q.dt[1]*(u128)BI::base;
			BI R2{BI::mul_small(R,2)};
			int c=R2<Ds?0:(Ds<R2?1:(M&1));
			if(c)M++;
			if(M>=(u128)1<<112){
				if(M==(u128)1<<112)return qmk(neg,1,0);
				M>>=1;
				return qmk(neg,1,M&(((u128)1<<112)-1));
			}
			return qmk(neg,0,M);
		}
	}
	FPTP
	inline i128 intok(Bint<T>&t,bool un){
		bool f{0};char c=rdc();
		while(c!=EOF&&(c<'0'||c>'9')){f=c=='-';c=rdc();}
		if(c==EOF)return -1;
		static V<i128>v{};static Z n{};
		v.clear();n=0;
		c=rdg(v,n,c);
		if(un)ungetc(c,stdin);
		bld(v,n,t,f);
		return 1;
	}
	FPTP
	inline bool rdf(T&t,bool un){
		int c=rdc();
		while(c!=EOF&&c!='-'&&c!='+'&&c!='.'&&c!='i'&&c!='I'&&c!='n'&&c!='N'&&(c<'0'||c>'9'))c=rdc();
		if(c==EOF)return 0;
		bool f{0};
		if(c=='-'||c=='+'){f=c=='-';c=rdc();}
		if(c==EOF)return 0;
		if(c=='i'||c=='I'){
			c=rdc();
			if(c=='n'||c=='N'){
				c=rdc();
				if(c=='f'||c=='F'){
					c=rdc();
					if(un&&c!=EOF)ungetc(c,stdin);
					t=f?-std::numeric_limits<T>::infinity():std::numeric_limits<T>::infinity();
					return 1;
				}
			}
			if(un&&c!=EOF)ungetc(c,stdin);
			return 0;
		}
		if(c=='n'||c=='N'){
			c=rdc();
			if(c=='a'||c=='A'){
				c=rdc();
				if(c=='n'||c=='N'){
					c=rdc();
					if(un&&c!=EOF)ungetc(c,stdin);
					t=f?-std::numeric_limits<T>::quiet_NaN():std::numeric_limits<T>::quiet_NaN();
					return 1;
				}
			}
			if(un&&c!=EOF)ungetc(c,stdin);
			return 0;
		}
		string s{};
		int fp{0},ex{0};
		bool any{0};
		while(c!=EOF&&c>='0'&&c<='9'){s+=c;any=1;c=rdc();}
		if(c=='.'){
			c=rdc();
			while(c!=EOF&&c>='0'&&c<='9'){s+=c;fp++;any=1;c=rdc();}
		}
		if(c=='e'||c=='E'){
			c=rdc();
			bool en{0};
			if(c=='-'||c=='+'){en=c=='-';c=rdc();}
			while(c!=EOF&&c>='0'&&c<='9'){if(ex<100000)ex=ex*10+(c-'0');c=rdc();}
			if(en)ex=-ex;
		}
		bool frac{c=='/'};
		if(un&&!frac&&c!=EOF)ungetc(c,stdin);
		if(!any){t=0;return 0;}
		Z b{0};
		while(b<s.size()&&s[b]=='0')b++;
		s.erase(0,b);
		if(s.empty())s="0";
		int drop{0};
		if(s.size()>38){
			drop=s.size()-38;
			char nx{s[38]};
			bool tl{false};
			for(Z i=39;i<s.size();i++)if(s[i]!='0'){tl=true;break;}
			s.resize(38);
			if(nx>'5'||(nx=='5'&&(tl||((s[37]-'0')&1)))){
				int i{37};
				for(;i>=0;i--){if(s[i]!='9'){s[i]++;break;}s[i]='0';}
				if(i<0){s="1"+s.substr(0,37);drop++;}
			}
		}
		u128 u{};
		for(char d:s)u=u*10+(d-'0');
		int q{ex-fp+drop};
		if constexpr(is_same_v<T,__float128>)t=qdec(f,u,q);
		else{
			T v=u;
			if(q>0)v*=fpf::pw10<T>(q);
			else if(q<0)v/=fpf::pw10<T>(-q);
			t=f?-v:v;
		}
		if(frac){T den{};if(!rdf(den,un))return 0;t=t/den;}
		return 1;
	}
	FPTP
	inline bool rdi(T&t,bool un){
		if constexpr(fpt<T>)return rdf(t,un);
		bool f{0};char c=rdc();t=0;
		while(c!=EOF&&(c<'0'||c>'9'))f=c=='-',c=rdc();
		if(c==EOF)return 0;
		while(c!=EOF&&!(c<'0'||c>'9'))t*=10,t+=c-'0',c=rdc();
		if(un)ungetc(c,stdin);
		t=f?-t:t;
		return 1;
	}
}
using dfpio::rdc;
template<typename T=i128>
Bint<T>readBint(){Bint<T>t{};dfpio::intok(t,0);return t;}
FPTP
i128 isread(Bint<T>&t){return dfpio::intok(t,1);}
FPTP
i128 isread(T&t){return dfpio::rdi(t,1)?1:-1;}
FPTP
i128 isread(Frac<T>&t){
	char c=rdc();
	while(c!=EOF&&(c<'0'||c>'9')&&c!='-'&&c!='+'&&c!='.')c=rdc();
	if(c==EOF)return -1;
	bool neg{0};
	if(c=='-'||c=='+'){neg=c=='-';c=rdc();}
	if(c==EOF)return -1;
	if constexpr(is_same_v<T,BI>){
		static V<i128>v{},w{};
		static Z vn{},wn{};
		v.clear();w.clear();vn=0;wn=0;
		if(c!='.')c=dfpio::rdg(v,vn,c);
		Z ipn{vn};
		if(c=='.'){c=rdc();c=dfpio::rdg(v,vn,c);}
		Z fpn{vn-ipn};
		if(c=='/'){c=rdc();c=dfpio::rdg(w,wn,c);}
		if(c!=EOF)ungetc(c,stdin);
		if(!vn)return -1;
		if(!wn)dfpio::dig(w,wn,1);
		for(Z i=0;i<fpn;i++)dfpio::dig(w,wn,0);
		BI a{},b{};
		dfpio::bld(v,vn,a,0);
		dfpio::bld(w,wn,b,0);
		t=Frac<T>{neg?-a:a,b};
		return 1;
	}else{
		string s{},d{};
		Z ipn{0};
		if(c!='.'){c=dfpio::rds(s,c);ipn=s.size();}
		if(c=='.'){c=rdc();c=dfpio::rds(s,c);}
		if(c=='/'){c=rdc();c=dfpio::rds(d,c);}
		if(c!=EOF)ungetc(c,stdin);
		if(s.empty())return -1;
		Z fpn{s.size()-ipn};
		if(!ipn)s="0"+s;
		auto from_str=[](const string&x){
			if constexpr(fpt<T>)return fpf::num<T>(x);
			else return T{x};
		};
		T num{from_str(s)},den{1};
		for(Z i=0;i<fpn;i++)den*=10;
		if(!d.empty())den*=from_str(d);
		t=Frac<T>{neg?-num:num,den};
		return 1;
	}
}
FPTP
i128 isread(cx<T>&z){
	char c=rdc();
	while(c!=EOF&&c!='('&&c!='i'&&c!='I'&&c!='-'&&c!='+'&&c!='.'&&(c<'0'||c>'9'))c=rdc();
	if(c==EOF)return -1;
	if(c=='('){
		T a{},b{};
		if(isread(a)==-1)return -1;
		c=rdc();
		if(c==','){if(isread(b)==-1)return -1;c=rdc();}
		while(c!=EOF&&c!=')')c=rdc();
		z={a,b};
		return 1;
	}
	bool n{0};
	if(c=='-'||c=='+'){n=c=='-';c=rdc();}
	if(c=='i'||c=='I'){
		c=rdc();
		if(c!=EOF)ungetc(c,stdin);
		z={T{},n?-T(1):T(1)};
		return 1;
	}
	if(c!=EOF)ungetc(c,stdin);
	T a{};
	if(isread(a)==-1)return -1;
	if(n)a=-a;
	c=rdc();
	if(c=='i'||c=='I'){z={T{},a};return 1;}
	if(c==','){
		T b{};
		if(isread(b)==-1)return -1;
		c=rdc();
		if(c==')')c=rdc();
		if(c!=EOF)ungetc(c,stdin);
		z={a,b};
		return 1;
	}
	if(c!='+'&&c!='-'){
		if(c!=EOF)ungetc(c,stdin);
		z={a,T{}};
		return 1;
	}
	bool m{c=='-'};
	c=rdc();
	if(c=='i'||c=='I'){z={a,m?-T(1):T(1)};return 1;}
	if(c!=EOF)ungetc(c,stdin);
	T b{};
	if(isread(b)==-1)return -1;
	if(m)b=-b;
	c=rdc();
	if(c=='i'||c=='I')z={a,b};
	else{
		if(c!=EOF)ungetc(c,stdin);
		z={a,T{}};
	}
	return 1;
}
FPTP
i128 isread(Poly<T>&p){
	V<T>v{};
	char c=rdc();
	while(true){
		if(c==EOF||c=='\n')break;
		ungetc(c,stdin);
		T t{};
		if(isread(t)==-1)break;
		v.push_back(move(t));
		c=rdc();
	}
	if(v.empty())return -1;
	p.dt=move(v);
	return 1;
}
FPTP
void write(const Bint<T>&t){
	static V<char>b{};
	Z n{t.dt.size()},p{};
	if(b.size()<n*20+2)b.resize(n*20+2);
	if(t.sig)b[p++]='-';
	dfpio::puti(b.data(),p,t.dt.back(),0);
	for(Z i=n-1;i--;)dfpio::puti(b.data(),p,t.dt[i],18);
	dfpio::wr(b.data(),p);
}
FPTP
void write(const T&t){
	T x{t};bool g{x<0};char b[64];int n{64};
	do{int d=x%10;if(d<0)d=-d;b[--n]='0'+d;x/=10;}while(x);
	if(g)putchar('-');
	dfpio::wr(b+n,(Z)(64-n));
}
inline void write(const float&x){printf("%.9g",x);}
inline void write(const double&x){printf("%.17g",x);}
inline void write(const long double&x){printf("%.21Lg",x);}
inline void write(const __float128&x){dfpio::qwr(x);}
FPTP
void write(const cx<T>&z){write(z.real());T b{z.imag()};if(b==0)b=T{};if(b<0)write(b);else{putchar('+');write(b);}putchar('i');}
template<>
inline void write(const char&c){putchar(c);}
FPTP
void write(const Frac<T>&t){
	static V<char>b{};
	if(t.a<0){putchar('-');write(Frac<T>{-t.a,t.b});return;}
	T r{};
	if constexpr(fpt<T>){
		T q{fpf::trunc(t.a/t.b)};
		write(q);
		r=fpf::mod(t.a,t.b);
	}else{
		write(t.a/t.b);
		r=t.a%t.b;
	}
	if(!r)return;
	putchar('.');
	std::map<T,Z>pos;
	V<int>dig;
	Z lim{fpt<T>?1024:(Z)-1};
	while(r&&!pos.count(r)&&dig.size()<lim){
		pos[r]=dig.size();
		r*=10;
		int d{};
		if constexpr(is_same_v<T,BI>)d=(r/t.b).dt[0];
		else{d=(r/t.b);if(d<0)d=0;if(d>9)d=9;}
		dig.push_back(d);
		if constexpr(fpt<T>)r=fpf::mod(r,t.b);else r%=t.b;
	}
	Z st{0},P{0},L{dig.size()};
	if(r&&pos.count(r)){st=pos[r];P=dig.size()-st;L=st+2*P;}
	if(b.size()<L)b.resize(L);
	for(Z i=0;i<L;i++)b[i]='0'+dig[i<dig.size()?i:st+(i-st)%P];
	dfpio::wr(b.data(),L);
	if(r)putchar('.'),putchar('.'),putchar('.');
}
FPTP
void write(const Poly<T>&p){
	putchar('P');putchar('{');
	for(Z i=0;i<p.dt.size();i++){if(i)putchar(',');write(p.dt[i]);}
	putchar('}');
}
FPTP
T read(){T t{};if constexpr(cxt<T>::value)isread(t);else dfpio::rdi(t,0);return t;}
template<typename T1,typename T2>
void write(const std::pair<T1,T2>&t){putchar('(');write(t.first);putchar(',');write(t.second);putchar(')');}
FPTP
void write(const V<T>&v){
	putchar('[');
	for(Z i=0;i<v.size();i++){if(i)putchar(',');write(v[i]);}
	putchar(']');
}
template<>
inline void write(const string&s){if(s.size())dfpio::wr(s.data(),s.size());}
FPTP
void writen(const T&f){write(f);putchar('\n');}
template<typename T,typename...Args>
void writen(const T&f,const Args&...args){write(f);putchar(' ');writen(args...);}
FPTP
void writeln(const T&f){write(f);putchar('\n');}
template<typename T,typename...Args>
void writeln(const T&f,const Args&...args){write(f);writeln(args...);}
inline i128 isread(string&s){
	s.clear();char c=rdc();
	while(c!=EOF&&(c==' '||c=='\t'||c=='\r'||c=='\n'))c=rdc();
	if(c==EOF)return -1;
	while(c!=EOF&&!(c==' '||c=='\t'||c=='\r'||c=='\n')){s+=c;c=rdc();}
	ungetc(c,stdin);
	return 1;
}
template<typename T1,typename T2>
i128 isread(std::pair<T1,T2>&p){if(isread(p.first)==-1)return -1;return isread(p.second);}
FPTP
i128 isread(V<T>&v){
	V<T>t{};
	char c=rdc();
	while(c!=EOF&&(c==' '||c=='\t'||c=='\r'))c=rdc();
	bool k{0};
	if(c=='['){k=1;c=rdc();}
	while(true){
		while(c!=EOF&&c!='\n'&&((k&&c!=']'&&(c==' '||c=='\t'||c=='\r'||c==','||c=='('||c==')'))||(!k&&(c==' '||c=='\t'||c=='\r'))))c=rdc();
		if(c==EOF||c=='\n'||(k&&c==']'))break;
		ungetc(c,stdin);
		T x{};
		if(isread(x)==-1)break;
		t.push_back(move(x));
		c=rdc();
	}
	if(k&&c==']'){char n=rdc();if(n!='\n')ungetc(n,stdin);}
	if(!k&&t.empty())return -1;
	v=move(t);
	return 1;
}
inline i128 isreadn(){return 0;}
template<typename T,typename...Args>
i128 isreadn(T&f,Args&...args){auto t1=isread(f),t2=isreadn(args...);return(t1==-1||t2==-1)?-1:t1+t2;}
#undef FPTP
#endif