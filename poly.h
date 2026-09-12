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
	static constexpr T base=1000000000;
	std::vector<T>dt{};
	bool sig{};
	Bint():dt{0},sig{}{}
	Bint(const T&t):dt{},sig{}{
		sig=t<0;
		__int128 x{t<0?-(__int128)t:(__int128)t};
		do{dt.push_back((T)(x%base));x/=base;}while(x);
	}
	Bint(const std::string&s):dt{},sig{}{
		for(char c:s){
			if(c=='-'){sig=1;continue;}
			__int128 carry{c-'0'};
			for(size_t i=0;carry||i<dt.size();i++){
				if(i==dt.size())dt.push_back({});
				__int128 cur{dt[i]*10+carry};
				dt[i]=(T)(cur%base);
				carry=cur/base;
			}
		}
		reload();
	}
	void reload(){
		while(dt.size()>1&&!dt.back())dt.pop_back();
		if(dt.empty())dt.push_back({});
		if(dt.size()==1&&!dt[0])sig=0;
	}
	explicit operator bool()const{return dt.size()>1||dt[0];}
	static bool abs_less(const Bint&x,const Bint&y){
		if(x.dt.size()!=y.dt.size())return x.dt.size()<y.dt.size();
		for(size_t i=x.dt.size();i--;)
			if(x.dt[i]!=y.dt[i])return x.dt[i]<y.dt[i];
		return false;
	}
	static Bint abs_add(const Bint&a,const Bint&b){
		Bint t{};
		t.dt.resize(std::max(a.dt.size(),b.dt.size())+1);
		__int128 carry{};
		for(size_t i=0;i<t.dt.size();i++){
			__int128 cur{carry};
			if(i<a.dt.size())cur+=a.dt[i];
			if(i<b.dt.size())cur+=b.dt[i];
			t.dt[i]=(T)(cur%base);
			carry=cur/base;
		}
		t.reload();
		return t;
	}
	static Bint abs_sub(const Bint&a,const Bint&b){
		Bint t{};
		t.dt.resize(a.dt.size());
		__int128 carry{};
		for(size_t i=0;i<a.dt.size();i++){
			__int128 cur{a.dt[i]-carry-(i<b.dt.size()?b.dt[i]:0)};
			if(cur<0)cur+=base,carry=1;
			else carry=0;
			t.dt[i]=(T)cur;
		}
		t.reload();
		return t;
	}
	static Bint mul_small(const Bint&a,__int128 c){
		Bint t{};
		t.dt.assign(a.dt.size()+1,{});
		__int128 carry{};
		for(size_t i=0;i<a.dt.size();i++){
			__int128 cur{a.dt[i]*c+carry};
			t.dt[i]=(T)(cur%base);
			carry=cur/base;
		}
		t.dt[a.dt.size()]=(T)carry;
		t.reload();
		return t;
	}
	static void divmod(const Bint&a,const Bint&b,Bint&q,Bint&r){
		q.dt.clear();r.dt.clear();q.sig=0;r.sig=0;
		for(size_t i=a.dt.size();i--;){
			r.dt.insert(r.dt.begin(),a.dt[i]);
			r.reload();
			__int128 c{};
			if(!abs_less(r,b)){
				__int128 R{r.dt.back()},B{b.dt.back()};
				if(r.dt.size()>1)R=R*base+r.dt[r.dt.size()-2];
				if(b.dt.size()>1)B=B*base+b.dt[b.dt.size()-2];
				__int128 lo{R/(B+1)},hi{base-1};
				if(lo<1)lo=1;
				while(lo<hi){
					__int128 mid{(lo+hi+1)>>1};
					if(abs_less(r,mul_small(b,mid)))hi=mid-1;
					else lo=mid;
				}
				c=lo;
				r=abs_sub(r,mul_small(b,c));
			}
			q.dt.insert(q.dt.begin(),(T)c);
		}
		q.reload();r.reload();
		if(q)q.sig=a.sig!=b.sig;
		if(r)r.sig=a.sig;
	}
	friend Bint operator+(const Bint&a,const Bint&b){
		Bint t{};
		if(a.sig==b.sig)t=abs_add(a,b),t.sig=a.sig;
		else if(abs_less(a,b))t=abs_sub(b,a),t.sig=b.sig;
		else t=abs_sub(a,b),t.sig=a.sig;
		t.reload();
		return t;
	}
	Bint&operator+=(const Bint&x){
		*this=*this+x;
		return *this;
	}
	friend Bint operator-(const Bint&a){
		Bint t{a};
		if(t)t.sig=!t.sig;
		return t;
	}
	friend Bint operator-(const Bint&a,const Bint&b){
		return a+(-b);
	}
	Bint&operator-=(const Bint&x){
		*this=*this-x;
		return *this;
	}
	friend Bint operator*(const Bint&a,const Bint&b){
		Bint t{};
		t.dt.assign(a.dt.size()+b.dt.size(),{});
		for(size_t i=0;i<a.dt.size();i++){
			__int128 carry{};
			for(size_t j=0;j<b.dt.size();j++){
				__int128 cur{(__int128)a.dt[i]*b.dt[j]+t.dt[i+j]+carry};
				t.dt[i+j]=(T)(cur%base);
				carry=cur/base;
			}
			for(size_t j=i+b.dt.size();carry;j++){
				__int128 cur{t.dt[j]+carry};
				t.dt[j]=(T)(cur%base);
				carry=cur/base;
			}
		}
		t.sig=a.sig!=b.sig;
		t.reload();
		return t;
	}
	Bint&operator*=(const Bint&x){
		*this=*this*x;
		return *this;
	}
	friend Bint operator/(const Bint&a,const Bint&b){
		Bint q{},r{};
		divmod(a,b,q,r);
		return q;
	}
	Bint&operator/=(const Bint&x){
		*this=*this/x;
		return *this;
	}
	friend Bint operator%(const Bint&a,const Bint&b){
		Bint q{},r{};
		divmod(a,b,q,r);
		return r;
	}
	Bint&operator%=(const Bint&x){
		*this=*this%x;
		return *this;
	}
	friend bool operator<(const Bint&a,const Bint&b){
		if(a.sig!=b.sig)return a.sig;
		return a.sig?abs_less(b,a):abs_less(a,b);
	}
	friend bool operator<=(const Bint&a,const Bint&b){
		return !(b<a);
	}
	friend bool operator==(const Bint&a,const Bint&b){
		return a.sig==b.sig&&a.dt==b.dt;
	}
	friend bool operator!=(const Bint&a,const Bint&b){
		return !(a==b);
	}
	friend bool operator>=(const Bint&a,const Bint&b){
		return !(a<b);
	}
	friend bool operator>(const Bint&a,const Bint&b){
		return b<a;
	}
};
template<typename T=Bint<__int128>>
class Frac{
	public:
	T a,b;
	Frac():a{},b{1}{}
	Frac(const T&_a,const T&_b=1):a{_a},b{_b}{simplify();}
	static constexpr long double eps=1e-12l;
	static T gcd(T x,T y){
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
			T g{gcd(a<0?-a:a,b)};
			if(g>1)a/=g,b/=g;
		}
	}
	explicit operator bool()const{return !!a;}
	friend Frac operator+(const Frac&x,const Frac&y){
		return {x.a*y.b+y.a*x.b,x.b*y.b};
	}
	Frac&operator+=(const Frac&x){
		*this=*this+x;
		return *this;
	}
	friend Frac operator-(const Frac&x){
		return {-x.a,x.b};
	}
	friend Frac operator-(const Frac&x,const Frac&y){
		return {x.a*y.b-y.a*x.b,x.b*y.b};
	}
	Frac&operator-=(const Frac&x){
		*this=*this-x;
		return *this;
	}
	friend Frac operator*(const Frac&x,const Frac&y){
		return {x.a*y.a,x.b*y.b};
	}
	Frac&operator*=(const Frac&x){
		*this=*this*x;
		return *this;
	}
	friend Frac operator/(const Frac&x,const Frac&y){
		return {x.a*y.b,x.b*y.a};
	}
	Frac&operator/=(const Frac&x){
		*this=*this/x;
		return *this;
	}
	friend bool operator<(const Frac&x,const Frac&y){
		return x.a*y.b<y.a*x.b;
	}
	friend bool operator<=(const Frac&x,const Frac&y){
		return !(y<x);
	}
	friend bool operator==(const Frac&x,const Frac&y){
		return eq(x.a*y.b,y.a*x.b);
	}
	friend bool operator!=(const Frac&x,const Frac&y){
		return !(x==y);
	}
	friend bool operator>=(const Frac&x,const Frac&y){
		return !(x<y);
	}
	friend bool operator>(const Frac&x,const Frac&y){
		return y<x;
	}
};
/**
 * 维护多项式 $\sum_{i=0}^{\texttt{dt.size()}}\texttt{dt[i]}x^i$ 相关的运算。
 * 有 operator+ operator+= operator-（一元运算符和二元运算符） operator-= operator* operator*= operator/ operator/= operator% operator%= operator< operator<= operator== operator!= operator>= operator>
 */
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
	static bool zero(const T&x){
		if constexpr(std::is_same_v<T,long double>)return std::fabs(x)<eps;
		return !x;
	}
	static bool eq(const T&a,const T&b){
		if constexpr(std::is_same_v<T,long double>)return std::fabs(a-b)<eps;
		return a==b;
	}
	void reload(){
		while(!dt.empty()&&zero(dt.back()))dt.pop_back();
		if(dt.empty())dt.push_back({});
	}
	T operator()(T x)const{
		T r{};
		for(size_t i=dt.size();i--;)r=r*x+dt[i];
		return r;
	}
	friend Poly operator+(const Poly&a,const Poly&b){
		Poly t{};
		t.dt.resize(std::max(a.dt.size(),b.dt.size()));
		for(size_t i=0;i<t.dt.size();i++){
			if(i<a.dt.size())t.dt[i]+=a.dt[i];
			if(i<b.dt.size())t.dt[i]+=b.dt[i];
		}
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
		Poly t{};
		t.dt.resize(a.dt.size());
		for(size_t i=0;i<a.dt.size();i++)t.dt[i]=-a.dt[i];
		t.reload();
		return t;
	}
	friend Poly operator-(const Poly&a,const Poly&b){
		Poly t{a};
		t-=b;
		return t;
	}
	Poly&operator-=(const Poly&b){
		if(dt.size()<b.dt.size())dt.resize(b.dt.size());
		for(size_t i=0;i<b.dt.size();i++)dt[i]-=b.dt[i];
		reload();
		return *this;
	}
	friend Poly operator*(const Poly&a,const Poly&b){
		Poly t{};
		t.dt.resize(a.dt.size()+b.dt.size()-1);
		for(size_t i=0;i<a.dt.size();i++)
			for(size_t j=0;j<b.dt.size();j++)
				t.dt[i+j]+=a.dt[i]*b.dt[j];
		t.reload();
		return t;
	}
	Poly&operator*=(const Poly&b){
		*this=*this*b;
		return *this;
	}
	friend Poly operator/(const Poly&a,const Poly&b){
		Poly t{a};
		t/=b;
		return t;
	}
	Poly&operator/=(const Poly&b){
		Poly t{},r{*this};
		t.dt.resize(r.dt.size());
		while(r.dt.size()>=b.dt.size()&&(r.dt.size()>1||!zero(r.dt[0]))){
			size_t sh=r.dt.size()-b.dt.size();
			t.dt[sh]=r.dt.back()/b.dt.back();
			for(size_t i=0;i<b.dt.size();i++)r.dt[sh+i]-=t.dt[sh]*b.dt[i];
			r.reload();
		}
		*this=std::move(t);
		reload();
		return *this;
	}
	friend Poly operator%(const Poly&a,const Poly&b){
		Poly t{a};
		t%=b;
		return t;
	}
	Poly&operator%=(const Poly&b){
		*this-=*this/b*b;
		reload();
		return *this;
	}
	friend bool operator<(const Poly&a,const Poly&b){
		Poly x{a},y{b};
		x.reload();y.reload();
		if(x.dt.size()!=y.dt.size())return x.dt.size()<y.dt.size();
		for(size_t i=x.dt.size();i--;)
			if(!eq(x.dt[i],y.dt[i]))return x.dt[i]<y.dt[i];
		return false;
	}
	friend bool operator<=(const Poly&a,const Poly&b){
		return !(b<a);
	}
	friend bool operator==(const Poly&a,const Poly&b){
		Poly x{a},y{b};
		x.reload();y.reload();
		if(x.dt.size()!=y.dt.size())return false;
		for(size_t i=0;i<x.dt.size();i++)
			if(!eq(x.dt[i],y.dt[i]))return false;
		return true;
	}
	friend bool operator!=(const Poly&a,const Poly&b){
		return !(a==b);
	}
	friend bool operator>=(const Poly&a,const Poly&b){
		return !(a<b);
	}
	friend bool operator>(const Poly&a,const Poly&b){
		return b<a;
	}
};
namespace std{
	template<typename T>
	Poly<T>gcd(const Poly<T>&a,const Poly<T>&b){
		Poly<T>x{a},y{b};
		x.reload();y.reload();
		while(y.dt.size()>1||!Poly<T>::zero(y.dt[0])){
			Poly<T>r{x%y};
			x=y;
			y=r;
		}
		if(x.dt.size()>1||!Poly<T>::zero(x.dt[0]))x/=Poly<T>{x.dt.back()};
		return x;
	}
}
template<typename T>
inline Poly<T>Deri(Poly<T>p){
	if(p.dt.size()==1){p.dt[0]={};return p;}
	for(size_t i=1;i<p.dt.size();i++)p.dt[i-1]=p.dt[i]*T(i);
	p.dt.pop_back();
	p.reload();
	return p;
}
template<typename T>
inline Poly<T>inter(Poly<T>p){
	p.dt.push_back({});
	for(size_t i=p.dt.size()-1;i--;)p.dt[i+1]=p.dt[i]/T(i+1);
	p.dt[0]={};
	p.reload();
	return p;
}
template<typename T=__int128>
inline Bint<T>readBint(){
	std::string s{};
	char c{getchar()};
	while(c!=EOF&&(c<'0'||c>'9')&&c!='-')c=getchar();
	if(c=='-'){s+='-';c=getchar();}
	while(c!=EOF&&!(c<'0'||c>'9')){s+=c;c=getchar();}
	return Bint<T>{s};
}
template<typename T>
inline __int128 isread(Bint<T>&t){
	bool f=0;char c=getchar();
	while(c!=EOF&&(c<'0'||c>'9')){f=c=='-';c=getchar();}
	if(c==EOF)return -1;
	std::string s{};
	while(c!=EOF&&!(c<'0'||c>'9')){s+=c;c=getchar();}
	ungetc(c,stdin);
	if(s.empty())return -1;
	t=Bint<T>{(f?"-":"")+s};
	return 1;
}
template<typename T>
inline __int128 isread(Frac<T>&t){
	auto from_str=[](const std::string&s){
		if constexpr(std::is_same_v<T,long double>)return std::stold(s);
		else return T{s};
	};
	char c{(char)getchar()};
	while(c!=EOF&&(c<'0'||c>'9')&&c!='-'&&c!='.')c=getchar();
	if(c==EOF)return -1;
	bool neg=0;
	if(c=='-'){
		neg=1;
		c=getchar();
	}
	if(c==EOF)return -1;
	std::string ip{},fp{},dp{};
	if(c!='.'){
		while(c!=EOF&&!(c<'0'||c>'9')){ip+=c;c=getchar();}
	}
	if(c=='.'){
		c=getchar();
		while(c!=EOF&&!(c<'0'||c>'9')){fp+=c;c=getchar();}
	}
	if(c=='/'){
		c=getchar();
		while(c!=EOF&&!(c<'0'||c>'9')){dp+=c;c=getchar();}
	}
	if(c!=EOF)ungetc(c,stdin);
	if(ip.empty()&&fp.empty())return -1;
	if(ip.empty())ip="0";
	T num{from_str(ip+fp)},den{1};
	for(size_t i=0;i<fp.size();i++)den*=10;
	if(!dp.empty())den*=from_str(dp);
	t=Frac<T>{neg?-num:num,den};
	return 1;
}
template<typename T>
inline __int128 isread(Poly<T>&p){
	std::vector<T>v{};
	while(true){
		char c=getchar();
		while(c!=EOF&&c!='\n'&&(c<'0'||c>'9')&&c!='-'&&c!='.')c=getchar();
		if(c==EOF||c=='\n')break;
		ungetc(c,stdin);
		T t{};
		if(isread(t)==-1)break;
		v.push_back(std::move(t));
	}
	if(v.empty())return -1;
	p.dt=std::move(v);
	return 1;
}
template<typename T>
inline void write(const Bint<T>&t){
	if(t.sig)putchar('-');
	printf("%d",(int)t.dt.back());
	for(size_t i=t.dt.size()-1;i--;)printf("%09d",(int)t.dt[i]);
}
template<typename T>
inline void write(const Frac<T>&t){
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
	if(!r){
		for(int c:dig)putchar('0'+c);
	}else{
		size_t st{pos[r]},L{dig.size()-st};
		for(size_t i=0;i<st+2*L;i++)
			putchar('0'+dig[i<dig.size()?i:st+(i-st)%L]);
		putchar('.');putchar('.');putchar('.');
	}
}
template<typename T>
inline void write(const Poly<T>&p){
	putchar('P');putchar('{');
	for(size_t i=0;i<p.dt.size();i++){
		if(i)putchar(',');
		write(p.dt[i]);
	}
	putchar('}');
}