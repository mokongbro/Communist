#ifndef VEC_POLY_HPP
#define VEC_POLY_HPP
#include<vector>
#include<map>
#include<string>
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<type_traits>
#include<numeric>
#include<initializer_list>
template<typename T=__int128>
class Bint{
	public:
	static constexpr T base=1000000000000000000;
	std::vector<T>dt{};
	bool sig{};
	Bint():dt{0},sig{}{}
	Bint(const T&t):dt{},sig{t<0}{
		__int128 x{t<0?-(__int128)t:(__int128)t};
		put(*this,(unsigned __int128)x);
	}
	Bint(const std::string&s):dt{},sig{}{
		size_t n{s.size()},d{};
		for(char c:s){if(c=='-'){sig=1;continue;}d++;}
		dt.assign((d+17)/18,{});
		size_t k{};
		unsigned long long w{},p{1};
		for(size_t i{n};i--;){
			char c{s[i]};
			if(c=='-')continue;
			w+=(unsigned long long)(c-'0')*p;
			p*=10;
			if(p==1000000000000000000ull){dt[k++]=(T)w;w=0;p=1;}
		}
		if(p!=1)dt[k++]=(T)w;
		reload();
	}
	void reload(){
		while(dt.size()>1&&!dt.back())dt.pop_back();
		if(dt.empty())dt.push_back({});
		if(dt.size()==1&&!dt[0])sig=0;
	}
	explicit operator bool()const{return dt.size()>1||dt[0];}
	static unsigned long long div9(unsigned __int128&x){
		unsigned long long r{},lo{(unsigned long long)x},hi{(unsigned long long)(x>>64)};
		unsigned long long w[4]{lo&4294967295ull,lo>>32,hi&4294967295ull,hi>>32},q[4]{};
		for(int i{3};i>=0;i--){unsigned long long c{r*4294967296ull+w[i]};q[i]=c/1000000000ull;r=c-q[i]*1000000000ull;}
		x=((unsigned __int128)(q[3]*4294967296ull+q[2])<<64)|(q[1]*4294967296ull+q[0]);
		return r;
	}
	static unsigned long long divbase(unsigned __int128&x){
		unsigned long long d0{div9(x)};
		unsigned long long d1{div9(x)};
		return d0+d1*1000000000ull;
	}
	static void put(Bint&t,unsigned __int128 x){
		t.dt.clear();
		if(!(x>>64)){unsigned long long y{(unsigned long long)x};do{t.dt.push_back((T)(y%1000000000000000000ull));y/=1000000000000000000ull;}while(y);}
		else do{t.dt.push_back((T)divbase(x));}while(x);
	}
	static Bint from128(__int128 v){Bint t{};t.sig=v<0;put(t,v<0?0-(unsigned __int128)v:(unsigned __int128)v);return t;}
	static bool smallu(const Bint&x,unsigned long long&v){
		if(x.dt.size()!=1)return 0;
		v=(unsigned long long)x.dt[0];
		return 1;
	}
	static bool smallv(const Bint&x,__int128&v){
		unsigned long long u{};
		if(!smallu(x,u))return 0;
		v=x.sig?-(__int128)u:(__int128)u;
		return 1;
	}
	static bool abs_less(const Bint&x,const Bint&y){
		if(x.dt.size()!=y.dt.size())return x.dt.size()<y.dt.size();
		for(size_t i=x.dt.size();i--;)if(x.dt[i]!=y.dt[i])return x.dt[i]<y.dt[i];
		return false;
	}
	static Bint abs_add(const Bint&a,const Bint&b){
		const std::vector<T>&x{a.dt.size()>=b.dt.size()?a.dt:b.dt};
		const std::vector<T>&y{a.dt.size()>=b.dt.size()?b.dt:a.dt};
		Bint t{};
		size_t n{x.size()},m{y.size()};
		t.dt.resize(n+1);
		unsigned long long c{};
		for(size_t i{0};i<m;i++){unsigned long long s{c+(unsigned long long)x[i]+(unsigned long long)y[i]};c=s/1000000000000000000ull;t.dt[i]=(T)(s-c*1000000000000000000ull);}
		for(size_t i{m};i<n;i++){unsigned long long s{c+(unsigned long long)x[i]};c=s/1000000000000000000ull;t.dt[i]=(T)(s-c*1000000000000000000ull);}
		t.dt[n]=(T)c;
		t.reload();
		return t;
	}
	static Bint abs_sub(const Bint&a,const Bint&b){
		Bint t{};
		size_t n{a.dt.size()},m{b.dt.size()};
		t.dt.resize(n);
		unsigned long long c{};
		size_t i{};
		for(;i<m&&i<n;i++){long long s{(long long)(unsigned long long)a.dt[i]-(long long)(unsigned long long)b.dt[i]-(long long)c};if(s<0)s+=1000000000000000000ll,c=1;else c=0;t.dt[i]=(T)s;}
		for(;i<n;i++){long long s{(long long)(unsigned long long)a.dt[i]-(long long)c};if(s<0)s+=1000000000000000000ll,c=1;else c=0;t.dt[i]=(T)s;}
		t.reload();
		return t;
	}
	static Bint mul_small(const Bint&a,__int128 c){Bint t{};t.dt.resize(a.dt.size()+1);__int128 k{};for(size_t i{0};i<a.dt.size();i++){__int128 v{(__int128)a.dt[i]*c+k};t.dt[i]=(T)(v%base);k=v/base;}t.dt[a.dt.size()]=(T)k;t.reload();return t;}
	static void cv(const T*a,size_t n,const T*b,size_t m,unsigned __int128*c){
		for(size_t i{0};i<n;i++){unsigned __int128 x{(unsigned long long)a[i]};for(size_t j{0};j<m;j++)c[i+j]+=x*(unsigned long long)b[j];}
	}
	static void nz(const unsigned __int128*c,size_t n,T*r){
		unsigned __int128 carry{};
		for(size_t i{0};i<n;i++){unsigned __int128 s{c[i]+carry};r[i]=(T)divbase(s);carry=s;}
	}
	static std::vector<T> sm(const T*a,size_t n,const T*b,size_t m){
		size_t L{n>m?n:m};
		std::vector<T>r(L+1);
		unsigned long long c{};
		for(size_t i{0};i<L;i++){
			unsigned long long s{c};
			if(i<n)s+=(unsigned long long)a[i];
			if(i<m)s+=(unsigned long long)b[i];
			c=s/1000000000000000000ull;
			r[i]=(T)(s-c*1000000000000000000ull);
		}
		r[L]=(T)c;
		while(r.size()>1&&!r.back())r.pop_back();
		return r;
	}
	static void sb(std::vector<T>&x,const std::vector<T>&y){
		__int128 c{};
		for(size_t i{0};i<x.size();i++){
			__int128 s{(__int128)x[i]-(i<y.size()?(__int128)y[i]:0)-c};
			if(s<0)s+=1000000000000000000,c=1;else c=0;
			x[i]=(T)s;
		}
		while(x.size()>1&&!x.back())x.pop_back();
	}
	static void as(std::vector<T>&r,const std::vector<T>&z,size_t k){
		if(r.size()<z.size()+k+1)r.resize(z.size()+k+1);
		unsigned long long c{};
		size_t i{0};
		for(;i<z.size();i++){unsigned long long s{(unsigned long long)r[k+i]+(unsigned long long)z[i]+c};c=s/1000000000000000000ull;r[k+i]=(T)(s-c*1000000000000000000ull);}
		for(size_t j{k+i};c;j++){unsigned long long s{(unsigned long long)r[j]+c};c=s/1000000000000000000ull;r[j]=(T)(s-c*1000000000000000000ull);}
	}
	static std::vector<T> mv(const T*a,size_t n,const T*b,size_t m){
		if(n<m){const T*p{a};a=b;b=p;size_t q{n};n=m;m=q;}
		std::vector<T>r{};
		if(m<=320){
			std::vector<unsigned __int128>c(n+m,0);
			cv(a,n,b,m,c.data());
			r.resize(n+m);
			nz(c.data(),n+m,r.data());
			while(r.size()>1&&!r.back())r.pop_back();
			return r;
		}
		size_t k{m>>1};
		std::vector<T>z0{mv(a,k,b,k)},z2{mv(a+k,n-k,b+k,m-k)};
		std::vector<T>s{sm(a,k,a+k,n-k)},t{sm(b,k,b+k,m-k)};
		std::vector<T>z1{mv(s.data(),s.size(),t.data(),t.size())};
		sb(z1,z0);
		sb(z1,z2);
		r=z0;
		as(r,z1,k);
		as(r,z2,2*k);
		while(r.size()>1&&!r.back())r.pop_back();
		return r;
	}
	static constexpr unsigned long long q1{998244353},q2{1004535809},q3{469762049},q12{1002772198720536577};
	static constexpr unsigned long long n1{998244351},n2{1004535807},n3{469762047},u1{932051910},u2{542374313},u3{460175152},w12{669690699},w3{354521948};
	static constexpr size_t th{2800};
	static unsigned long long mm(unsigned long long a,unsigned long long b,unsigned long long p,unsigned long long np){
		unsigned long long t{a*b},m{(t&4294967295ull)*np&4294967295ull},r{(t+m*p)>>32};
		return r>=p?r-p:r;
	}
	static unsigned long long mp(unsigned long long a,unsigned long long e,unsigned long long p,unsigned long long np,unsigned long long r2){
		unsigned long long t{mm(a,r2,p,np)},r{mm(1,r2,p,np)};
		for(;e;e>>=1){if(e&1)r=mm(r,t,p,np);t=mm(t,t,p,np);}
		return r;
	}
	static void rt(std::vector<unsigned long long>&t,size_t n,unsigned long long p,unsigned long long np,unsigned long long r2,bool iv){
		t.resize(n);
		for(size_t l{2};l<=n;l<<=1){
			unsigned long long e{iv?p-1-(p-1)/l:(p-1)/l},wl{mp(3,e,p,np,r2)},w{mm(1,r2,p,np)};
			for(size_t j{0};j<l/2;j++){t[l/2+j]=w;w=mm(w,wl,p,np);}
		}
	}
	static void tf(std::vector<unsigned long long>&f,size_t n,const std::vector<unsigned long long>&t,unsigned long long p,unsigned long long np,unsigned long long r2,bool iv){
		for(size_t i{1},j{};i<n;i++){
			size_t b{n>>1};
			for(;j&b;b>>=1)j^=b;
			j^=b;
			if(i<j){unsigned long long q{f[i]};f[i]=f[j];f[j]=q;}
		}
		for(size_t l{2};l<=n;l<<=1)for(size_t i{0};i<n;i+=l)for(size_t j{0};j<l/2;j++){
			unsigned long long u{f[i+j]},v{mm(f[i+j+l/2],t[l/2+j],p,np)},s{u+v};
			f[i+j]=s>=p?s-p:s;
			f[i+j+l/2]=u>=v?u-v:u+p-v;
		}
		if(iv){
			unsigned long long s{mp(n%p,p-2,p,np,r2)};
			for(size_t i{0};i<n;i++)f[i]=mm(f[i],s,p,np);
		}
	}
	static void nv(std::vector<unsigned long long>&f,std::vector<unsigned long long>&g,const unsigned long long*a,size_t n,const unsigned long long*b,size_t m,size_t N,unsigned long long p,unsigned long long np,unsigned long long r2,std::vector<unsigned long long>&o){
		f.resize(N);g.resize(N);
		for(size_t i{0};i<n;i++){unsigned long long v{a[i]};while(v>=p)v-=p;f[i]=mm(v,r2,p,np);}
		for(size_t i{n};i<N;i++)f[i]=0;
		for(size_t i{0};i<m;i++){unsigned long long v{b[i]};while(v>=p)v-=p;g[i]=mm(v,r2,p,np);}
		for(size_t i{m};i<N;i++)g[i]=0;
		std::vector<unsigned long long>t{};
		rt(t,N,p,np,r2,0);
		tf(f,N,t,p,np,r2,0);tf(g,N,t,p,np,r2,0);
		for(size_t i{0};i<N;i++)f[i]=mm(f[i],g[i],p,np);
		rt(t,N,p,np,r2,1);
		tf(f,N,t,p,np,r2,1);
		o.resize(n+m-1);
		for(size_t i{0};i<o.size();i++)o[i]=mm(f[i],1,p,np);
	}
	static std::vector<T> nm(const Bint&a,const Bint&b){
		std::vector<unsigned long long>ha(a.dt.size()*2),hb(b.dt.size()*2),f,g,r1,r2,r3;
		size_t n{half(a,ha.data())},m{half(b,hb.data())},L{n+m-1},N{1};
		while(N<L)N<<=1;
		nv(f,g,ha.data(),n,hb.data(),m,N,q1,n1,u1,r1);
		nv(f,g,ha.data(),n,hb.data(),m,N,q2,n2,u2,r2);
		std::vector<unsigned long long>c(L);
		for(size_t i{0};i<L;i++)c[i]=(unsigned long long)((unsigned __int128)r1[i]+q1*((unsigned __int128)(r2[i]>=r1[i]?r2[i]-r1[i]:r2[i]+q2-r1[i])*w12%q2));
		nv(f,g,ha.data(),n,hb.data(),m,N,q3,n3,u3,r3);
		std::vector<unsigned long long>d{};
		d.reserve(L+4);
		unsigned __int128 cy{};
		for(size_t i{0};i<L;i++){
			unsigned long long v{(unsigned long long)((unsigned __int128)(r3[i]>=c[i]%q3?r3[i]-c[i]%q3:r3[i]+q3-c[i]%q3)*w3%q3)};
			unsigned __int128 s{(unsigned __int128)c[i]+(unsigned __int128)q12*v+cy};
			d.push_back(div9(s));
			cy=s;
		}
		while(cy)d.push_back(div9(cy));
		std::vector<T>r((d.size()+1)/2);
		for(size_t i{0};i<r.size();i++)r[i]=(T)(d[2*i]+(2*i+1<d.size()?d[2*i+1]*1000000000ull:0));
		while(r.size()>1&&!r.back())r.pop_back();
		return r;
	}
	friend Bint operator*(const Bint&a,const Bint&b){
		__int128 x{},y{};
		if(smallv(a,x)&&smallv(b,y))return from128(x*y);
		size_t n{a.dt.size()},m{b.dt.size()};
		Bint t{};
		t.sig=a.sig!=b.sig;
		if(n<=32&&m<=32){
			unsigned __int128 c[64];
			for(size_t i{0};i<n+m;i++)c[i]=0;
			cv(a.dt.data(),n,b.dt.data(),m,c);
			t.dt.resize(n+m);
			nz(c,n+m,t.dt.data());
		}else if((n<m?n:m)>=th&&n+m<=1048576)t.dt=nm(a,b);
		else t.dt=mv(a.dt.data(),n,b.dt.data(),m);
		t.reload();
		return t;
	}
	static size_t half(const Bint&x,unsigned long long*h){
		size_t k{};
		for(size_t i{0};i<x.dt.size();i++){
			unsigned long long v{(unsigned long long)x.dt[i]};
			h[k++]=v%1000000000ull;
			h[k++]=v/1000000000ull;
		}
		while(k>1&&!h[k-1])k--;
		return k;
	}
	static constexpr size_t bth=32;
	static unsigned long long nf(unsigned long long t){return t<500000000000000000ull?(500000000000000000ull+(t-1))/t:1;}
	static void tb(Bint&x){while(x.dt.size()>1&&!x.dt.back())x.dt.pop_back();if(x.dt.empty())x.dt.push_back({});}
	static Bint lo(const Bint&x,size_t k){
		Bint t{};
		size_t n{x.dt.size()<k?x.dt.size():k};
		if(n)t.dt.assign(x.dt.begin(),x.dt.begin()+n);
		return t;
	}
	static Bint hi(const Bint&x,size_t k){
		Bint t{};
		if(x.dt.size()>k)t.dt.assign(x.dt.begin()+k,x.dt.end());
		return t;
	}
	static Bint cc(const Bint&h,const Bint&l,size_t k){
		Bint t{};
		t.dt.assign(k,{});
		size_t n{l.dt.size()<k?l.dt.size():k};
		for(size_t i{0};i<n;i++)t.dt[i]=l.dt[i];
		for(size_t i{0};i<h.dt.size();i++)t.dt.push_back(h.dt[i]);
		tb(t);
		return t;
	}
	static Bint bl(const Bint&x,size_t o,size_t k){
		Bint t{};
		size_t b{o<x.dt.size()?o:x.dt.size()},e{o+k<x.dt.size()?o+k:x.dt.size()};
		if(b<e)t.dt.assign(x.dt.begin()+b,x.dt.begin()+e);
		return t;
	}
	static Bint ds(const Bint&x,unsigned long long d){
		Bint t{};
		t.dt.resize(x.dt.size());
		unsigned long long c{};
		for(size_t i{x.dt.size()};i--;){
			unsigned __int128 v{(unsigned __int128)c*1000000000000000000ull+(unsigned long long)x.dt[i]};
			unsigned long long q{(unsigned long long)(v/d)};
			t.dt[i]=(T)q;
			c=(unsigned long long)(v-(unsigned __int128)q*d);
		}
		tb(t);
		return t;
	}
	static void kd(const Bint&A,const Bint&C,Bint&q,Bint&r){
		size_t n{A.dt.size()},m{C.dt.size()};
		if(m==1){
			unsigned long long d{(unsigned long long)C.dt[0]},k{};
			q.dt.resize(n);
			for(size_t i{n};i--;){
				unsigned __int128 x{(unsigned __int128)k*1000000000000000000ull+(unsigned long long)A.dt[i]};
				unsigned long long t{(unsigned long long)(x/d)};
				q.dt[i]=(T)t;
				k=(unsigned long long)(x-(unsigned __int128)t*d);
			}
			tb(q);
			r.dt.assign(1,(T)k);
			return;
		}
		unsigned long long d{nf((unsigned long long)C.dt[m-1])},k{};
		std::vector<unsigned long long>u(n+1),v(m);
		if(d==1){
			for(size_t i{0};i<n;i++)u[i]=(unsigned long long)A.dt[i];
			for(size_t i{0};i<m;i++)v[i]=(unsigned long long)C.dt[i];
		}else{
			for(size_t i{0};i<n;i++){unsigned __int128 p{(unsigned __int128)(unsigned long long)A.dt[i]*d+k};u[i]=(unsigned long long)(p%1000000000000000000ull);k=(unsigned long long)(p/1000000000000000000ull);}
			u[n]=k;k=0;
			for(size_t i{0};i<m;i++){unsigned __int128 p{(unsigned __int128)(unsigned long long)C.dt[i]*d+k};v[i]=(unsigned long long)(p%1000000000000000000ull);k=(unsigned long long)(p/1000000000000000000ull);}
		}
		q.dt.assign(n-m+1,(T)0);
		for(size_t j{n-m+1};j--;){
			unsigned __int128 nm{(unsigned __int128)u[j+m]*1000000000000000000ull+u[j+m-1]};
			unsigned long long qh{(unsigned long long)(nm/(unsigned __int128)v[m-1])},rh{(unsigned long long)(nm%(unsigned __int128)v[m-1])};
			if(qh>=1000000000000000000ull){qh=999999999999999999ull;rh=(unsigned long long)(nm-(unsigned __int128)qh*v[m-1]);}
			while(rh<1000000000000000000ull&&(unsigned __int128)qh*v[m-2]>(unsigned __int128)rh*1000000000000000000ull+u[j+m-2])qh--,rh+=v[m-1];
			unsigned long long c2{},br{};
			for(size_t i{0};i<m;i++){
				unsigned __int128 p{(unsigned __int128)qh*v[i]+c2};
				c2=(unsigned long long)(p/1000000000000000000ull);
				unsigned long long pl{(unsigned long long)(p-c2*1000000000000000000ull)},t{u[j+i]-pl-br};
				br=t>>63;
				u[j+i]=t+br*1000000000000000000ull;
			}
			unsigned long long t0{u[j+m]-c2-br};
			br=t0>>63;
			u[j+m]=t0+br*1000000000000000000ull;
			while(br){
				qh--;
				unsigned long long g{};
				for(size_t i{0};i<m;i++){unsigned long long s{u[j+i]+v[i]+g};if(s>=1000000000000000000ull){s-=1000000000000000000ull;g=1;}else g=0;u[j+i]=s;}
				u[j+m]+=g;
				if(u[j+m]>=1000000000000000000ull)u[j+m]-=1000000000000000000ull;
				br=u[j+m]!=0;
			}
			q.dt[j]=(T)qh;
		}
		tb(q);
		r.dt.resize(m);
		if(d==1){for(size_t i{0};i<m;i++)r.dt[i]=(T)u[i];}
		else{unsigned long long g{};for(size_t i{m};i--;){unsigned __int128 x{(unsigned __int128)g*1000000000000000000ull+u[i]};unsigned long long t{(unsigned long long)(x/d)};r.dt[i]=(T)t;g=(unsigned long long)(x-(unsigned __int128)t*d);}}
		tb(r);
	}
	static void bz(const Bint&A,const Bint&C,Bint&q,Bint&r){
		size_t n{A.dt.size()},m{C.dt.size()};
		if(n<m){q=Bint{};r=A;return;}
		if(m<=bth){kd(A,C,q,r);return;}
		size_t h{m>>1},K{(n-m+1+h-1)/h};
		Bint R{bl(A,h*(K-1),n-h*(K-1))};
		q.dt.assign(h*K,(T)0);
		for(size_t i{K};i--;){
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
			for(size_t j{0};j<qa.dt.size()&&h*i+j<q.dt.size();j++)q.dt[h*i+j]=qa.dt[j];
			R=ra;
		}
		tb(q);
		r=R;
	}
	static void divmod(const Bint&a,const Bint&b,Bint&q,Bint&r){
		Bint A{a},C{b};
		A.sig=0;C.sig=0;
		A.reload();C.reload();
		if(!(bool)A||abs_less(A,C)){
			q=Bint{};
			r=A;
			if(r)r.sig=a.sig;
			return;
		}
		size_t m{C.dt.size()};
		if(m<=bth)kd(A,C,q,r);
		else{
			unsigned long long d{nf((unsigned long long)C.dt[m-1])};
			if(d!=1)A=mul_small(A,(__int128)d),C=mul_small(C,(__int128)d);
			bz(A,C,q,r);
			if(d!=1)r=ds(r,d);
		}
		q.reload();r.reload();
		if(q)q.sig=a.sig!=b.sig;
		if(r)r.sig=a.sig;
	}
	friend Bint operator+(const Bint&a,const Bint&b){
		__int128 x{},y{};
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
		__int128 x{},y{};
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
template<typename U>struct fracb_isbint{static constexpr bool value=false;};
template<typename U>struct fracb_isbint<Bint<U>>{static constexpr bool value=true;};
template<typename T=Bint<__int128>>
class Frac{
	public:
	T a,b;
	Frac():a{},b{1}{}
	Frac(const T&_a,const T&_b=1):a{_a},b{_b}{simplify();}
	template<typename U=T,typename=std::enable_if_t<fracb_isbint<U>::value>>
	Frac(const __int128&x):a{x},b{1}{}
	static constexpr long double eps=1e-12l;
	static bool one(const T&x){
		if constexpr(fracb_isbint<T>::value)return !x.sig&&x.dt.size()==1&&x.dt[0]==1;
		return x==1;
	}
	static unsigned __int128 wd(const T&x){return x.dt.size()>1?(unsigned __int128)x.dt[1]*(unsigned __int128)T::base+(unsigned __int128)x.dt[0]:(unsigned __int128)x.dt[0];}
	static T mkw(unsigned __int128 v){
		T t{};
		t.dt.clear();
		do{t.dt.push_back(v%(unsigned __int128)T::base);v/=(unsigned __int128)T::base;}while(v);
		return t;
	}
	static unsigned __int128 gw(unsigned __int128 x,unsigned __int128 y){
		while(y){unsigned __int128 t{x%y};x=y;y=t;}
		return x;
	}
	static T ed(const T&a,unsigned __int128 m){
		T t{a};unsigned __int128 c{0};
		for(size_t i=t.dt.size();i--;){unsigned __int128 v{c*(unsigned __int128)T::base+(unsigned __int128)t.dt[i]};t.dt[i]=(__int128)(v/m);c=v%m;}
		while(t.dt.size()>1&&!t.dt.back())t.dt.pop_back();
		return t;
	}
	static T dv(const T&a,const T&g){
		if(a.dt==g.dt)return a.sig?T{-1}:T{1};
		if(g.dt.size()==1)return ed(a,(unsigned __int128)g.dt[0]);
		if(a.dt.size()<=2&&g.dt.size()<=2){T q{mkw(wd(a)/wd(g))};q.sig=a.sig;return q;}
		return a/g;
	}
	static void hl(T&x){__int128 c{0};
		for(size_t i=x.dt.size();i--;){__int128 v{(__int128)x.dt[i]+c*T::base};x.dt[i]=v/2;c=v&1;}
		while(x.dt.size()>1&&!x.dt.back())x.dt.pop_back();
	}
	static void dl(T&x){
		__int128 c{0};
		for(size_t i=0;i<x.dt.size();i++){__int128 v{(__int128)x.dt[i]*2+c};x.dt[i]=v%T::base;c=v/T::base;}
		if(c)x.dt.push_back(c);
	}
	static void sb(T&x,const T&y){__int128 c{0};
		for(size_t i=0;i<x.dt.size();i++){__int128 v{(__int128)x.dt[i]-c-(i<y.dt.size()?(__int128)y.dt[i]:0)};if(v<0)v+=T::base,c=1;else c=0;x.dt[i]=v;}
		while(x.dt.size()>1&&!x.dt.back())x.dt.pop_back();
	}
	static T bg2(T x,T y){
		size_t z{0},w{0};
		while(!(x.dt[0]&1)){hl(x);z++;}
		while(!(y.dt[0]&1)){hl(y);w++;}
		if(z>w)z=w;
		for(;;){if(T::abs_less(x,y))std::swap(x,y);if(x.dt==y.dt)break;sb(x,y);while(!(x.dt[0]&1))hl(x);}
		while(z--)dl(x);return x;
	}
	static T gcd(T x,T y){
		if constexpr(fracb_isbint<T>::value){
			for(;;){
				if(T::abs_less(x,y))std::swap(x,y);
				if(!y)return x;
				if(y.dt.size()==1){
					unsigned __int128 m{(unsigned __int128)y.dt[0]},r{0};if(m<2)return mkw(1);
					for(size_t i=x.dt.size();i--;)r=(r*(unsigned __int128)T::base+(unsigned __int128)x.dt[i])%m;
					return mkw(gw(r,m));
				}
				if(x.dt.size()<=2)return mkw(gw(wd(x),wd(y)));
				if(x.dt.size()<=4)return bg2(x,y);
				T r{x%y};x=y;y=r;
			}
		}
		if(x<0)x=-x;
		if(y<0)y=-y;
		while(y){T t{x%y};x=y;y=t;}
		return x;
	}
	static bool eq(const T&x,const T&y){
		if constexpr(std::is_same_v<T,long double>)return std::fabs(x-y)<eps;
		return x==y;
	}
	void simplify(){
		if(b<0)a=-a,b=-b;
		if constexpr(std::is_same_v<T,long double>)return;
		else{
			if(one(b))return;
			T g{gcd(a,b)};
			if(g>1)a=dv(a,g),b=dv(b,g);
		}
	}
	explicit operator bool()const{return !!a;}
	static Frac qt(const T&p,const T&q,const T&r,const T&s){
		T g{gcd(p,s)};
		if(g>1){T h{gcd(r,q)};if(h>1)return {dv(p,g)*dv(r,h),dv(q,h)*dv(s,g)};return {dv(p,g)*r,q*dv(s,g)};}
		T h{gcd(r,q)};
		if(h>1)return {p*dv(r,h),dv(q,h)*s};
		return {p*r,q*s};
	}
	static Frac ad(const Frac&x,const Frac&y,bool m){
		if constexpr(fracb_isbint<T>::value){
			if(x.b.dt.size()+y.b.dt.size()>=3&&x.b&&y.b){
				T g{gcd(x.b,y.b)};
				if(g>1){T g1{dv(x.b,g)},g2{dv(y.b,g)};return {m?x.a*g2-y.a*g1:x.a*g2+y.a*g1,g1*y.b};}
			}
		}
		return {m?x.a*y.b-y.a*x.b:x.a*y.b+y.a*x.b,x.b*y.b};
	}
	friend Frac operator+(const Frac&x,const Frac&y){return ad(x,y,0);}
	friend Frac operator-(const Frac&x,const Frac&y){return ad(x,y,1);}
	friend Frac operator-(const Frac&x){return {-x.a,x.b};}
	friend Frac operator*(const Frac&x,const Frac&y){
		if constexpr(fracb_isbint<T>::value)if(x.a.dt.size()+y.a.dt.size()>=5)return qt(x.a,x.b,y.a,y.b);
		return {x.a*y.a,x.b*y.b};
	}
	friend Frac operator/(const Frac&x,const Frac&y){
		if constexpr(fracb_isbint<T>::value)if(x.a.dt.size()+y.b.dt.size()>=5)return qt(x.a,x.b,y.b,y.a);
		return {x.a*y.b,x.b*y.a};
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
/**
 * 维护多项式 $\sum_{i=0}^{\texttt{dt.size()}}\texttt{dt[i]}x^i$ 相关的运算。
 * 有 operator+ operator+= operator-（一元运算符和二元运算符） operator-= operator* operator*= operator/ operator/= operator% operator%= operator< operator<= operator== operator!= operator>= operator>
 */
template<typename T>struct pex{using type=T;static constexpr bool v=false;};
template<typename U>struct pex<Frac<Bint<U>>>{using type=Bint<U>;static constexpr bool v=true;};
template<typename T=Frac<Bint<__int128>>>
class Poly{
	public:
	std::vector<T>dt;
	void assign(size_t _n,const T&x){dt.assign(_n,x);}
	T&operator[](size_t i){return dt[i];}
	const T&operator[](size_t i)const{return dt[i];}
	auto begin(){return dt.begin();}
	auto end(){return dt.end();}
	auto begin()const{return dt.begin();}
	auto end()const{return dt.end();}
	Poly(){}
	Poly(const std::vector<T>&_dt):dt{_dt}{}
	Poly(std::initializer_list<T>_dt):dt{_dt}{}
	Poly(const T&t):dt{t}{}
	static constexpr long double eps=1e-12l;
	static bool zero(const T&x){if constexpr(std::is_same_v<T,long double>)return std::fabs(x)<eps;return !x;}
	static bool eq(const T&a,const T&b){if constexpr(std::is_same_v<T,long double>)return std::fabs(a-b)<eps;return a==b;}
	void reload(){while(!dt.empty()&&zero(dt.back()))dt.pop_back();if(dt.empty())dt.push_back({});}
	T operator()(T x)const{T r{};for(size_t i=dt.size();i--;)r=r*x+dt[i];return r;}
	template<typename A>
	static A ig(A x,A y){if(x.sig)x=-x;if(y.sig)y=-y;while(y){A t{x%y};x=y;y=t;}return x;}
	template<typename A>
	static A lc(const std::vector<T>&v){A l{1};for(const T&x:v){if(x.b==l||(x.b.dt.size()==1&&x.b.dt[0]==1))continue;l=l/ig(l,x.b)*x.b;if(l.dt.size()>8)return A{};}return l;}
	template<typename A>
	static std::vector<A> conv(const std::vector<A>&a,const std::vector<A>&b){
		std::vector<A>u{a},v{b},c;
		if(u.size()<v.size())u.swap(v);
		c.resize(u.size()+v.size()-1);
		if(v.size()<=32){for(size_t i=0;i<u.size();i++)for(size_t j=0;j<v.size();j++)c[i+j]+=u[i]*v[j];return c;}
		size_t h{v.size()/2};
		std::vector<A>u1(u.begin()+h,u.end()),v1(v.begin()+h,v.end()),z0,z1,z2;
		u.resize(h);v.resize(h);
		z0=conv(u,v);z2=conv(u1,v1);
		u.resize(u1.size());v.resize(v1.size());
		for(size_t i=0;i<u1.size();i++)u[i]+=u1[i];for(size_t i=0;i<v1.size();i++)v[i]+=v1[i];
		z1=conv(u,v);
		for(size_t i=0;i<z0.size();i++)z1[i]-=z0[i];for(size_t i=0;i<z2.size();i++)z1[i]-=z2[i];
		for(size_t i=0;i<z1.size();i++){if(i<z0.size())c[i]+=z0[i];c[i+h]+=z1[i];c[i+2*h]+=z2[i];}
		return c;
	}
	static bool smul(const Poly&a,const Poly&b,Poly&t){
		using A=typename pex<T>::type;
		A la{lc<A>(a.dt)},lb{lc<A>(b.dt)};
		if(!la||!lb)return false;
		std::vector<A>ua(a.dt.size()),ub(b.dt.size());
		for(size_t i=0;i<ua.size();i++)ua[i]=a.dt[i].b==la?a.dt[i].a:a.dt[i].a*(la/a.dt[i].b);
		for(size_t i=0;i<ub.size();i++)ub[i]=b.dt[i].b==lb?b.dt[i].a:b.dt[i].a*(lb/b.dt[i].b);
		std::vector<A>uc{conv(ua,ub)};
		A den{la*lb};
		t.dt.resize(uc.size());
		if(den.dt.size()==1&&den.dt[0]==1&&!den.sig)for(size_t i=0;i<uc.size();i++)t.dt[i].a=uc[i];
		else for(size_t i=0;i<uc.size();i++)t.dt[i]=T{uc[i],den};
		t.reload();
		return true;
	}
	static void pdiv(const std::vector<T>&a,const std::vector<T>&b,std::vector<T>&q){
		size_t m{b.size()},sz{a.size()};
		std::vector<T>r{a};
		q.assign(sz,T{});
		while(sz>=m&&(sz>1||!zero(r[0]))){
			size_t sh{sz-m};
			T c{r[sz-1]/b[m-1]};
			q[sh]=c;
			if constexpr(pex<T>::v){
				for(size_t i=0;i+1<m;i++){
					T&ri{r[sh+i]};
					auto cb{c.b*b[i].b};
					ri=T{ri.a*cb-c.a*b[i].a*ri.b,ri.b*cb};
				}
				r[sh+m-1]=T{};
			}else for(size_t i=0;i<m;i++)r[sh+i]-=c*b[i];
			while(sz>1&&zero(r[sz-1]))sz--;
		}
	}
	friend Poly operator+(const Poly&a,const Poly&b){
		Poly t{};
		t.dt.resize(std::max(a.dt.size(),b.dt.size()));
		for(size_t i=0;i<t.dt.size();i++){if(i<a.dt.size())t.dt[i]+=a.dt[i];if(i<b.dt.size())t.dt[i]+=b.dt[i];}
		t.reload();
		return t;
	}
	Poly&operator+=(const Poly&b){
		if(dt.size()<b.dt.size())dt.resize(b.dt.size());
		for(size_t i=0;i<b.dt.size();i++)dt[i]+=b.dt[i];
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
		for(size_t i=0;i<b.dt.size();i++)dt[i]-=b.dt[i];
		reload();
		return *this;
	}
	friend Poly operator*(const Poly&a,const Poly&b){
		Poly t{};
		if constexpr(pex<T>::v)if(!a.dt.empty()&&!b.dt.empty()&&smul(a,b,t))return t;
		t.dt.resize(a.dt.size()+b.dt.size()-1);
		for(size_t i=0;i<a.dt.size();i++)for(size_t j=0;j<b.dt.size();j++)t.dt[i+j]+=a.dt[i]*b.dt[j];
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
		for(size_t i=x.dt.size();i--;)if(!eq(x.dt[i],y.dt[i]))return x.dt[i]<y.dt[i];
		return false;
	}
	friend bool operator<=(const Poly&a,const Poly&b){return !(b<a);}
	friend bool operator==(const Poly&a,const Poly&b){
		Poly x{a},y{b};
		x.reload();y.reload();
		if(x.dt.size()!=y.dt.size())return false;
		for(size_t i=0;i<x.dt.size();i++)if(!eq(x.dt[i],y.dt[i]))return false;
		return true;
	}
	friend bool operator!=(const Poly&a,const Poly&b){return !(a==b);}
	friend bool operator>=(const Poly&a,const Poly&b){return !(a<b);}
	friend bool operator>(const Poly&a,const Poly&b){return b<a;}
};
namespace fpg{
using u64=unsigned long long;
using FB=Bint<__int128>;
using FF=Frac<FB>;
using FP=Poly<FF>;
static const unsigned pf[6]={1000000007u,1000000009u,998244353u,1004535809u,469762049u,167772161u};
inline u64 pw(u64 a,u64 e,u64 p){u64 r{1};a%=p;while(e){if(e&1)r=r*a%p;a=a*a%p;e>>=1;}return r;}
inline unsigned mb(const FB&x,unsigned p){
	u64 r{},b=(u64)(FB::base%(__int128)p);
	for(size_t i=x.dt.size();i--;)r=(r*b+(u64)(x.dt[i]%(__int128)p))%p;
	if(x.sig&&r)r=p-r;
	return (unsigned)r;
}
inline bool red(const FP&a,unsigned p,std::vector<unsigned>&v){
	v.resize(a.dt.size());
	for(size_t i=0;i<a.dt.size();i++){
		unsigned n=mb(a.dt[i].a,p),d=mb(a.dt[i].b,p);
		if(!d)return 0;
		v[i]=(unsigned)((u64)n*pw(d,p-2,p)%p);
	}
	while(v.size()>1&&!v.back())v.pop_back();
	return v.size()>1||v[0];
}
inline void gcdp(std::vector<unsigned>a,std::vector<unsigned>b,unsigned p,std::vector<unsigned>&g){
	while(!(b.size()==1&&!b[0])){
		std::vector<unsigned>r{a};
		unsigned v=(unsigned)pw(b.back(),p-2,p);
		while(r.size()>=b.size()&&!(r.size()==1&&!r[0])){
			size_t sh=r.size()-b.size();
			unsigned c=(unsigned)((u64)r.back()*v%p);
			if(c)for(size_t i=0;i<b.size();i++){
				unsigned t=(unsigned)((u64)c*b[i]%p);
				r[sh+i]=r[sh+i]>=t?(unsigned)(r[sh+i]-t):(unsigned)(r[sh+i]+p-t);
			}
			while(r.size()>1&&!r.back())r.pop_back();
		}
		a=b;b=r;
	}
	unsigned v=(unsigned)pw(a.back(),p-2,p);
	for(auto&x:a)x=(unsigned)((u64)x*v%p);
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
	std::vector<FB>X{};
	int k{-1};
	for(int it=0;it<6;it++){
		unsigned p=pf[it];
		if(!mb(la,p)||!mb(lb,p))continue;
		std::vector<unsigned>A,B,g;
		if(!red(a,p,A)||!red(b,p,B))continue;
		if(A.size()!=a.dt.size()||B.size()!=b.dt.size())continue;
		gcdp(A,B,p,g);
		int kk=(int)g.size()-1;
		if(!kk){o=FP{1};return 1;}
		if(k<0)k=kk;
		if(kk!=k)continue;
		unsigned c=mb(la,p);
		std::vector<unsigned>r(g.size());
		for(size_t j=0;j<g.size();j++)r[j]=(unsigned)((u64)c*g[j]%p);
		if(X.empty()){
			X.resize(g.size());
			for(size_t j=0;j<g.size();j++)X[j]=FB{(__int128)r[j]};
			M=FB{(__int128)p};
		}else{
			unsigned iv=(unsigned)pw(mb(M,p),p-2,p);
			for(size_t j=0;j<g.size();j++){
				unsigned x=mb(X[j],p),t=x<=r[j]?(unsigned)(r[j]-x):(unsigned)(r[j]+p-x);
				unsigned u=(unsigned)((u64)t*iv%p);
				if(u)X[j]+=M*FB{(__int128)u};
			}
			M*=FB{(__int128)p};
		}
		if(it){
			FP q{};
			q.dt.resize(X.size());
			for(size_t j=0;j<X.size();j++){
				FB v{X[j]};
				if(v*FB{(__int128)2}>M)v-=M;
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
	template<typename T>
	Poly<T>gcd(const Poly<T>&a,const Poly<T>&b){
		if constexpr(std::is_same_v<T,Frac<Bint<__int128>>>){Poly<T>r{};if(fpg::fast(a,b,r))return r;}
		Poly<T>x{a},y{b};
		x.reload();y.reload();
		while(y.dt.size()>1||!Poly<T>::zero(y.dt[0])){Poly<T>r{x%y};x=std::move(y);y=std::move(r);}
		if(x.dt.size()>1||!Poly<T>::zero(x.dt[0]))x/=Poly<T>{x.dt.back()};
		return x;
	}
}
template<typename T>
inline Poly<T>Deri(Poly<T>p){
	if(p.dt.size()==1){p.dt[0]={};return p;}
	if constexpr(pex<T>::v){
		using A=typename pex<T>::type;
		for(size_t i=1;i<p.dt.size();i++)p.dt[i-1]=T{p.dt[i].a*A(i),p.dt[i].b};
	}else for(size_t i=1;i<p.dt.size();i++)p.dt[i-1]=p.dt[i]*T(i);
	p.dt.pop_back();
	p.reload();
	return p;
}
template<typename T>
inline Poly<T>inter(Poly<T>p){
	p.dt.push_back({});
	if constexpr(pex<T>::v){
		using A=typename pex<T>::type;
		for(size_t i=p.dt.size()-1;i--;)p.dt[i+1]=T{p.dt[i].a,p.dt[i].b*A(i+1)};
	}else for(size_t i=p.dt.size()-1;i--;)p.dt[i+1]=p.dt[i]/T(i+1);
	p.dt[0]={};
	p.reload();
	return p;
}
template<typename T>
inline Poly<T>radic(const Poly<T>&a){
	if(a.dt.size()==1&&Poly<T>::zero(a.dt[0]))return a;
	Poly<T>r{a/std::gcd(a,Deri(a))};
	if(r.dt.size()==1&&Poly<T>::zero(r.dt[0]))return r;
	return r/r.dt.back();
}
template<typename T1,typename T2>
inline T1 fpow(T1 a,T2 b){
	T1 ans{1};
	while(b){
		if(b&1)ans*=a;
		b>>=1;
		if(b)a*=a;
	}
	return ans;
}
template<typename T1,typename T2>
inline T1 fpow(T1 a,T2 b,T1 p){
	T1 ans{1};
	a%=p;
	while(b){
		if(b&1)(ans*=a)%=p;
		b>>=1;
		if(b)(a*=a)%=p;
	}
	return ans;
}
namespace dfpio{
#ifdef _getchar_nolock
	inline int rdc(){return _getchar_nolock();}
#else
	inline int rdc(){return getchar();}
#endif
	inline void dig(std::vector<__int128>&v,size_t&n,int d){
		if(n%18==0)v.push_back(d);
		else v.back()=v.back()*10+d;
		n++;
	}
	template<typename T>
	inline void bld(const std::vector<__int128>&v,size_t n,Bint<T>&t,bool neg){
		if(!n){t.dt.assign(1,{});t.sig=0;return;}
		size_t r{n%18};if(!r)r=18;
		__int128 cur{v.back()},p{1};
		for(size_t j=0;j<r;j++)p*=10;
		t.dt.clear();t.dt.reserve(v.size()+1);
		for(size_t i=v.size()-1;i--;){
			cur+=v[i]*p;
			t.dt.push_back((T)(cur%1000000000000000000));
			cur/=1000000000000000000;
		}
		t.dt.push_back((T)cur);
		t.sig=neg;
		t.reload();
	}
	inline char rdg(std::vector<__int128>&v,size_t&n,char c){
		while(c!=EOF&&!(c<'0'||c>'9')){dig(v,n,c-'0');c=(char)rdc();}
		return c;
	}
	inline char rds(std::string&s,char c){
		while(c!=EOF&&!(c<'0'||c>'9')){s+=c;c=(char)rdc();}
		return c;
	}
	inline void wr(const char*b,size_t n){fwrite(b,1,n,stdout);}
	inline void puti(char*b,size_t&p,unsigned long long m,int w){
		static const char t[]="00010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899";
		char d[20];int n{0};
		while(m>=100){unsigned r{(unsigned)(m%100)};m/=100;d[n++]=t[r*2+1];d[n++]=t[r*2];}
		if(m>9){d[n++]=t[m*2+1];d[n++]=t[m*2];}else d[n++]=(char)('0'+m);
		int q{w-n};if(q<0)q=0;
		while(q--)b[p++]='0';
		while(n)b[p++]=d[--n];
	}
	template<typename T>
	inline __int128 intok(Bint<T>&t,bool un){
		bool f{0};char c{(char)rdc()};
		while(c!=EOF&&(c<'0'||c>'9')){f=c=='-';c=(char)rdc();}
		if(c==EOF)return -1;
		static std::vector<__int128>v{};static size_t n{};
		v.clear();n=0;
		c=rdg(v,n,c);
		if(un)ungetc(c,stdin);
		bld(v,n,t,f);
		return 1;
	}
	template<typename T>
	inline bool rdi(T&t,bool un){
		bool f{0};char c{(char)rdc()};t=0;
		while(c!=EOF&&(c<'0'||c>'9'))f=c=='-',c=(char)rdc();
		if(c==EOF)return 0;
		while(c!=EOF&&!(c<'0'||c>'9'))t*=10,t+=c-'0',c=(char)rdc();
		if(un)ungetc(c,stdin);
		t=f?-t:t;
		return 1;
	}
}
template<typename T=__int128>
inline Bint<T>readBint(){Bint<T>t{};dfpio::intok(t,0);return t;}
template<typename T>
inline __int128 isread(Bint<T>&t){return dfpio::intok(t,1);}
template<typename T>
inline __int128 isread(Frac<T>&t){
	char c{(char)dfpio::rdc()};
	while(c!=EOF&&(c<'0'||c>'9')&&c!='-'&&c!='.')c=(char)dfpio::rdc();
	if(c==EOF)return -1;
	bool neg{0};
	if(c=='-'){neg=1;c=(char)dfpio::rdc();}
	if(c==EOF)return -1;
	if constexpr(std::is_same_v<T,Bint<__int128>>){
		static std::vector<__int128>v{},w{};
		static size_t vn{},wn{};
		v.clear();w.clear();vn=0;wn=0;
		if(c!='.')c=dfpio::rdg(v,vn,c);
		size_t ipn{vn};
		if(c=='.'){c=(char)dfpio::rdc();c=dfpio::rdg(v,vn,c);}
		size_t fpn{vn-ipn};
		if(c=='/'){c=(char)dfpio::rdc();c=dfpio::rdg(w,wn,c);}
		if(c!=EOF)ungetc(c,stdin);
		if(!vn)return -1;
		if(!wn)dfpio::dig(w,wn,1);
		for(size_t i=0;i<fpn;i++)dfpio::dig(w,wn,0);
		Bint<__int128>a{},b{};
		dfpio::bld(v,vn,a,0);
		dfpio::bld(w,wn,b,0);
		t=Frac<T>{neg?-a:a,b};
		return 1;
	}else{
		std::string s{},d{};
		size_t ipn{0};
		if(c!='.'){c=dfpio::rds(s,c);ipn=s.size();}
		if(c=='.'){c=(char)dfpio::rdc();c=dfpio::rds(s,c);}
		if(c=='/'){c=(char)dfpio::rdc();c=dfpio::rds(d,c);}
		if(c!=EOF)ungetc(c,stdin);
		if(s.empty())return -1;
		size_t fpn{s.size()-ipn};
		if(!ipn)s="0"+s;
		auto from_str=[](const std::string&x){
			if constexpr(std::is_same_v<T,long double>)return std::stold(x);
			else return T{x};
		};
		T num{from_str(s)},den{1};
		for(size_t i=0;i<fpn;i++)den*=10;
		if(!d.empty())den*=from_str(d);
		t=Frac<T>{neg?-num:num,den};
		return 1;
	}
}
template<typename T>
inline __int128 isread(Poly<T>&p){
	std::vector<T>v{};
	char c{(char)dfpio::rdc()};
	while(true){
		while(c!=EOF&&c!='\n'&&(c<'0'||c>'9')&&c!='-'&&c!='.')c=(char)dfpio::rdc();
		if(c==EOF||c=='\n')break;
		ungetc(c,stdin);
		T t{};
		if(isread(t)==-1)break;
		v.push_back(std::move(t));
		c=(char)dfpio::rdc();
	}
	if(v.empty())return -1;
	p.dt=std::move(v);
	return 1;
}
template<typename T>
inline void write(const Bint<T>&t){
	static std::vector<char>b{};
	size_t n{t.dt.size()},p{};
	if(b.size()<n*20+2)b.resize(n*20+2);
	if(t.sig)b[p++]='-';
	dfpio::puti(b.data(),p,(unsigned long long)t.dt.back(),0);
	for(size_t i=n-1;i--;)dfpio::puti(b.data(),p,(unsigned long long)t.dt[i],18);
	dfpio::wr(b.data(),p);
}
template<typename T>
inline void write(const Frac<T>&t){
	static std::vector<char>b{};
	if(t.a<0){putchar('-');write(Frac<T>{-t.a,t.b});return;}
	write(t.a/t.b);
	T r{t.a%t.b};
	if(!r)return;
	putchar('.');
	std::map<T,size_t>pos;
	std::vector<int>dig;
	while(r&&!pos.count(r)){
		pos[r]=dig.size();
		r*=10;
		int d{};
		if constexpr(std::is_same_v<T,Bint<__int128>>)d=(int)(r/t.b).dt[0];
		else d=(int)(r/t.b);
		dig.push_back(d);
		r%=t.b;
	}
	size_t st{0},P{0},L{dig.size()};
	if(r){st=pos[r];P=dig.size()-st;L=st+2*P;}
	if(b.size()<L)b.resize(L);
	for(size_t i=0;i<L;i++)b[i]=(char)('0'+dig[i<dig.size()?i:st+(i-st)%P]);
	dfpio::wr(b.data(),L);
	if(r)putchar('.'),putchar('.'),putchar('.');
}
template<typename T>
inline void write(const Poly<T>&p){
	putchar('P');putchar('{');
	for(size_t i=0;i<p.dt.size();i++){if(i)putchar(',');write(p.dt[i]);}
	putchar('}');
}
template<typename T>
inline T read(){T t{};dfpio::rdi(t,0);return t;}
template<typename T>
inline void write(const T&t){
	T x{t};bool g{x<0};char b[64];int n{64};
	do{int d{(int)(x%10)};if(d<0)d=-d;b[--n]=(char)('0'+d);x/=10;}while(x);
	if(g)putchar('-');
	dfpio::wr(b+n,(size_t)(64-n));
}
template<>
inline void write(const char&c){putchar(c);}
template<>
inline void write(const long double&x){printf("%.14Lf",x);}
template<typename T1,typename T2>
inline void write(const std::pair<T1,T2>&t){putchar('(');write(t.first);putchar(',');write(t.second);putchar(')');}
template<typename T>
inline void write(const std::vector<T>&v){
	putchar('[');
	for(size_t i=0;i<v.size();i++){if(i)putchar(',');write(v[i]);}
	putchar(']');
}
template<>
inline void write(const std::string&s){if(s.size())dfpio::wr(s.data(),s.size());}
template<typename T>
inline void writen(const T&f){write(f);putchar('\n');}
template<typename T,typename...Args>
inline void writen(const T&f,const Args&...args){write(f);putchar(' ');writen(args...);}
template<typename T>
inline void writeln(const T&f){write(f);putchar('\n');}
template<typename T,typename...Args>
inline void writeln(const T&f,const Args&...args){write(f);writeln(args...);}
template<typename T>
inline __int128 isread(T&t){return dfpio::rdi(t,1)?1:-1;}
inline __int128 isread(std::string&s){
	s.clear();char c{(char)dfpio::rdc()};
	while(c!=EOF&&(c==' '||c=='\t'||c=='\r'||c=='\n'))c=(char)dfpio::rdc();
	if(c==EOF)return -1;
	while(c!=EOF&&!(c==' '||c=='\t'||c=='\r'||c=='\n')){s+=c;c=(char)dfpio::rdc();}
	ungetc(c,stdin);
	return 1;
}
template<typename T1,typename T2>
inline __int128 isread(std::pair<T1,T2>&p){if(isread(p.first)==-1)return -1;return isread(p.second);}
template<typename T>
inline __int128 isread(std::vector<T>&v){
	std::vector<T>t{};
	char c{(char)dfpio::rdc()};
	while(c!=EOF&&(c==' '||c=='\t'||c=='\r'))c=(char)dfpio::rdc();
	bool k{0};
	if(c=='['){k=1;c=(char)dfpio::rdc();}
	while(true){
		while(c!=EOF&&c!='\n'&&((k&&c!=']'&&(c==' '||c=='\t'||c=='\r'||c==','||c=='('||c==')'))||(!k&&(c==' '||c=='\t'||c=='\r'))))c=(char)dfpio::rdc();
		if(c==EOF||c=='\n'||(k&&c==']'))break;
		ungetc(c,stdin);
		T x{};
		if(isread(x)==-1)break;
		t.push_back(std::move(x));
		c=(char)dfpio::rdc();
	}
	if(k&&c==']'){char n{(char)dfpio::rdc()};if(n!='\n')ungetc(n,stdin);}
	if(!k&&t.empty())return -1;
	v=std::move(t);
	return 1;
}
inline __int128 isreadn(){return 0;}
template<typename T,typename...Args>
inline __int128 isreadn(T&f,Args&...args){auto t1=isread(f),t2=isreadn(args...);return(t1==-1||t2==-1)?-1:t1+t2;}
#endif