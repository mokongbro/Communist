#ifndef D_SEG_H
#define D_SEG_H
using i128=__int128;
using u128=unsigned __int128;
#include<stdio.h>
template<typename T=__int128>
inline T read(){T t{};bool f=0;char c=getchar();
	while(c!=EOF&&(c<'0'||c>'9'))f=c=='-',c=getchar();
	while(c!=EOF&&!(c<'0'||c>'9'))t*=10,t+=c-'0',c=getchar();
	return f?-t:t;
}
template<typename T>
inline void write(T t){
	if(t<0)putchar('-'),t=-t;
	if(t>9)write(t/10);
	putchar(t%10+'0');
}
#include<vector>
#include<algorithm>
inline i128 mod(i128 x,i128 m){return x%m;}
const i128 INF=(static_cast<i128>(1)<<126)+((static_cast<i128>(1)<<126)-1);
const i128 NEG=-INF-1;
class Info{
	public:
	i128 l{},r{},d{};
	Info(const i128&_l=0,const i128&_r=0,const i128&_d=0):l(_l),r(_r),d(_d){};
};
class Seg{
	class nd{
		public:
		i128 sum{},tagp{},tagm{1},lc{NEG},rc{NEG};
		nd():sum{},tagp{},tagm{1},lc{NEG},rc{NEG}{}
	};
	std::vector<nd>dt;
	i128 L{},R{},MOD{INF};
	public:
	Seg():L{},R{},MOD{INF}{}
	void build(i128 n){L=1;R=n;if(!dt.size())dt.emplace_back();}
	void build(i128 _L,i128 _R){L=_L;R=_R;if(!dt.size())dt.emplace_back();}
	void setMod(i128 _m){MOD=_m;}
	void clear(){dt.clear();L=0;R=0;MOD=INF;}
	size_t size(){
		return dt.size();
	}
	const i128&newNd(i128&u,const i128&val=0){
		if(u!=NEG)return u;
		u=size();
		dt.emplace_back();
		dt.back().sum=val;
		return u;
	}
	void upd(const i128&u,const i128&s,const i128&t,const i128&m){
		auto lc=dt[u].lc,rc=dt[u].rc;
		dt[lc].sum=mod(dt[lc].sum*dt[u].tagm+(m-s+1)*dt[u].tagp,MOD);
		dt[rc].sum=mod(dt[rc].sum*dt[u].tagm+(t-m)*dt[u].tagp,MOD);
		dt[lc].tagm=mod(dt[lc].tagm*dt[u].tagm,MOD);
		dt[rc].tagm=mod(dt[rc].tagm*dt[u].tagm,MOD);
		dt[lc].tagp=mod(dt[lc].tagp*dt[u].tagm+dt[u].tagp,MOD);
		dt[rc].tagp=mod(dt[rc].tagp*dt[u].tagm+dt[u].tagp,MOD);
		dt[u].tagm=1;dt[u].tagp=0;
	}
	void add(const i128&d,const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return;
		if(l<=s&&t<=r){
			dt[u].sum=mod(dt[u].sum+(t-s+1)*d,MOD);
			dt[u].tagp=mod(dt[u].tagp+d,MOD);
			return;
		}
		auto m=(s+t)>>1;
		newNd(dt[u].lc,0);newNd(dt[u].rc,0);
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(l<=m)add(d,lc,l,r,s,m);
		if(r>m)add(d,rc,l,r,m+1,t);
		dt[u].sum=mod(dt[lc].sum+dt[rc].sum,MOD);
	}
	void add(const i128&l,const i128&r,const i128&x){
		add(x,0,l,r,L,R);
	}
	void mul(const i128&d,const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return;
		if(l<=s&&t<=r){
			dt[u].sum=mod(dt[u].sum*d,MOD);
			dt[u].tagp=mod(dt[u].tagp*d,MOD);
			dt[u].tagm=mod(dt[u].tagm*d,MOD);
			return;
		}
		auto m=(s+t)>>1;
		newNd(dt[u].lc);newNd(dt[u].rc);
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(l<=m)mul(d,lc,l,r,s,m);
		if(r>m)mul(d,rc,l,r,m+1,t);
		dt[u].sum=mod(dt[lc].sum+dt[rc].sum,MOD);
	}
	void mul(const i128&l,const i128&r,const i128&x){
		mul(x,0,l,r,L,R);
	}
	i128 query(const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return 0;
		if(l<=s&&t<=r)return dt[u].sum;
		i128 ans{};
		newNd(dt[u].lc);newNd(dt[u].rc);
		auto lc=dt[u].lc,rc=dt[u].rc;
		auto m=(s+t)>>1;
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		if(l<=m)ans=mod(ans+query(lc,l,r,s,m),MOD);
		if(r>m)ans=mod(ans+query(rc,l,r,m+1,t),MOD);
		return ans;
	}
	i128 query(const i128&l,const i128&r){
		return query(0,l,r,L,R);
	}
	//求第k大的数，此函数假设Seg是一个桶。
	i128 find(const i128&k,const i128&u,const i128&s,const i128&t){
		if(s==t)return s;
		newNd(dt[u].lc);newNd(dt[u].rc);
		auto lc=dt[u].lc,rc=dt[u].rc;
		auto m=(s+t)>>1;
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		if(k<=0)return s;
		if(k<=dt[rc].sum)return find(k,rc,m+1,t);
		if(k<=dt[rc].sum+dt[lc].sum)return find(k-dt[rc].sum,lc,s,m);
		return t;
	}
	i128 find(const i128&k){
		return find(k,0,L,R);
	}
	void reset(const i128&x,const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return;
		if(l<=s&&t<=r){
			dt[u].sum=mod((t-s+1)*x,MOD);
			dt[u].tagm=0;dt[u].tagp=mod(x,MOD);
			return;
		}
		auto m=(s+t)>>1;
		newNd(dt[u].lc);newNd(dt[u].rc);
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(l<=m)reset(x,lc,l,r,s,m);
		if(r>m)reset(x,rc,l,r,m+1,t);
		dt[u].sum=mod(dt[lc].sum+dt[rc].sum,MOD);
	}
	void reset(const i128&l,const i128&r,const i128&x){
		reset(x,0,l,r,L,R);
	}
	Seg&operator+=(Info p){
		add(p.l,p.r,p.d);
		return *this;
	}
	Seg&operator*=(Info p){
		mul(p.l,p.r,p.d);
		return *this;
	}
	i128 operator()(i128 l,i128 r){
		return query(l,r);
	}
	i128 operator[](const i128&k){
		return find(k);
	}
	void operator=(Info p){
		reset(p.l,p.r,p.d);
	}
};
class maxmin_Seg{
	class nd{
		public:
		i128 mx{NEG},mn{INF},tagp{},tagm{1},lc{NEG},rc{NEG};
		nd():mx{NEG},mn{INF},tagp{},tagm{1},lc{NEG},rc{NEG}{}
		void swapm(){if(mx<mn)std::swap(mx,mn);}
	};
	std::vector<nd>dt;
	i128 L{},R{},MOD{INF};
	public:
	maxmin_Seg():L{},R{},MOD{INF}{}
	void build(i128 n){L=1;R=n;if(!dt.size())dt.emplace_back();}
	void build(i128 _L,i128 _R){L=_L;R=_R;if(!dt.size())dt.emplace_back();}
	void setMod(i128 _m){MOD=_m;}
	void clear(){dt.clear();L=0;R=0;MOD=INF;}
	size_t size(){
		return dt.size();
	}
	const i128&newNd(i128&u,const i128&pmx=0,const i128&pmn=0){
		if(u!=NEG)return u;
		u=size();
		dt.emplace_back();
		dt.back().mx=pmx;dt.back().mn=pmn;
		return u;
	}
	void upd(const i128&u,const i128&s,const i128&t,const i128&m){
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(lc!=NEG){
			dt[lc].mx=mod(dt[lc].mx*dt[u].tagm+dt[u].tagp,MOD);
			dt[lc].mn=mod(dt[lc].mn*dt[u].tagm+dt[u].tagp,MOD);
			dt[lc].swapm();
			dt[lc].tagm=mod(dt[lc].tagm*dt[u].tagm,MOD);
			dt[lc].tagp=mod(dt[lc].tagp*dt[u].tagm+dt[u].tagp,MOD);
		}
		if(rc!=NEG){
			dt[rc].mx=mod(dt[rc].mx*dt[u].tagm+dt[u].tagp,MOD);
			dt[rc].mn=mod(dt[rc].mn*dt[u].tagm+dt[u].tagp,MOD);
			dt[rc].swapm();
			dt[rc].tagm=mod(dt[rc].tagm*dt[u].tagm,MOD);
			dt[rc].tagp=mod(dt[rc].tagp*dt[u].tagm+dt[u].tagp,MOD);
		}
		dt[u].tagm=1;dt[u].tagp=0;
	}
	void add(const i128&d,const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return;
		if(l<=s&&t<=r){
			dt[u].mx=mod(dt[u].mx+d,MOD);
			dt[u].mn=mod(dt[u].mn+d,MOD);
			dt[u].tagp=mod(dt[u].tagp+d,MOD);
			return;
		}
		auto m=(s+t)>>1;
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		newNd(dt[u].lc,dt[u].mx,dt[u].mn);newNd(dt[u].rc,dt[u].mx,dt[u].mn);
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(l<=m)add(d,lc,l,r,s,m);
		if(r>m)add(d,rc,l,r,m+1,t);
		dt[u].mx=std::max(dt[lc].mx,dt[rc].mx);
		dt[u].mn=std::min(dt[lc].mn,dt[rc].mn);
	}
	void add(const i128&l,const i128&r,const i128&x){
		add(x,0,l,r,L,R);
	}
	void mul(const i128&d,const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return;
		if(l<=s&&t<=r){
			dt[u].mx=mod(dt[u].mx*d,MOD);
			dt[u].mn=mod(dt[u].mn*d,MOD);
			dt[u].tagp=mod(dt[u].tagp*d,MOD);
			dt[u].tagm=mod(dt[u].tagm*d,MOD);
			dt[u].swapm();
			return;
		}
		auto m=(s+t)>>1;
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		newNd(dt[u].lc,dt[u].mx,dt[u].mn);newNd(dt[u].rc,dt[u].mx,dt[u].mn);
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(l<=m)mul(d,lc,l,r,s,m);
		if(r>m)mul(d,rc,l,r,m+1,t);
		dt[u].mx=std::max(dt[lc].mx,dt[rc].mx);
		dt[u].mn=std::min(dt[lc].mn,dt[rc].mn);
	}
	void mul(const i128&l,const i128&r,const i128&x){
		mul(x,0,l,r,L,R);
	}
	i128 querymx(const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return NEG;
		if(l<=s&&t<=r)return dt[u].mx;
		i128 ans=NEG;
		auto m=(s+t)>>1;
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		newNd(dt[u].lc,dt[u].mx,dt[u].mn);newNd(dt[u].rc,dt[u].mx,dt[u].mn);
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(l<=m)(ans=std::max(ans,querymx(lc,l,r,s,m)));
		if(r>m)(ans=std::max(ans,querymx(rc,l,r,m+1,t)));
		return ans;
	}
	i128 querymx(const i128&l,const i128&r){
		return querymx(0,l,r,L,R);
	}
	i128 querymn(const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return INF;
		if(l<=s&&t<=r)return dt[u].mn;
		i128 ans=INF;
		auto m=(s+t)>>1;
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		newNd(dt[u].lc,dt[u].mx,dt[u].mn);newNd(dt[u].rc,dt[u].mx,dt[u].mn);
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(l<=m)(ans=std::min(ans,querymn(lc,l,r,s,m)));
		if(r>m)(ans=std::min(ans,querymn(rc,l,r,m+1,t)));
		return ans;
	}
	i128 querymn(const i128&l,const i128&r){
		return querymn(0,l,r,L,R);
	}
	void reset(const i128&x,const i128&u,const i128&l,const i128&r,const i128&s,const i128&t){
		if(t<l||r<s)return;
		if(l<=s&&t<=r){
			dt[u].mx=mod(x,MOD);dt[u].mn=mod(x,MOD);
			dt[u].tagm=0;dt[u].tagp=mod(x,MOD);
			return;
		}
		auto m=(s+t)>>1;
		if(dt[u].tagm!=1||dt[u].tagp!=0)upd(u,s,t,m);
		newNd(dt[u].lc,dt[u].mx,dt[u].mn);newNd(dt[u].rc,dt[u].mx,dt[u].mn);
		auto lc=dt[u].lc,rc=dt[u].rc;
		if(l<=m)reset(x,lc,l,r,s,m);
		if(r>m)reset(x,rc,l,r,m+1,t);
		dt[u].mx=std::max(dt[lc].mx,dt[rc].mx);
		dt[u].mn=std::min(dt[lc].mn,dt[rc].mn);
	}
	void reset(const i128&l,const i128&r,const i128&x){
		reset(x,0,l,r,L,R);
	}
	maxmin_Seg&operator+=(Info p){
		add(p.l,p.r,p.d);
		return *this;
	}
	maxmin_Seg&operator*=(Info p){
		mul(p.l,p.r,p.d);
		return *this;
	}
	i128 operator()(i128 l,i128 r,i128 op){
		if(op>0)return querymx(l,r);
		else if(op<0)return querymn(l,r);
		else return 0;
	}
	void operator=(Info p){
		reset(p.l,p.r,p.d);
	}
};
#endif