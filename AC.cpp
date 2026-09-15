#include<iostream>
#include<cstdio>
#include<vector>
#include<string>
#include<array>
#include<queue>
#include<map>
template<typename T=__int128>
inline T read(){T t{};bool f=0;char c=getchar();
	while(c!=EOF&&(c<'0'||c>'9'))f=c=='-',c=getchar();
	while(c!=EOF&&!(c<'0'||c>'9'))t*=10,t+=c-'0',c=getchar();
	return f?-t:t;
}
template<typename T>
inline void write(const T&t){
	if(t<0)return putchar('-'),write(-t);
	if(t>9)write(t/10);
	putchar(t%10+'0');
}
template<>
inline void write(const char&c){
	putchar(c);
}
template<>
inline void write(const long double&x){
	printf("%.14Lf",x);
}
template<typename T1,typename T2>
inline void write(const std::pair<T1,T2>&t){
	putchar('(');write(t.first);putchar(',');write(t.second);putchar(')');
}
template<typename T>
inline void write(const std::vector<T>&v){
	putchar('[');
	for(std::size_t i=0;i<v.size();i++){
		if(i)putchar(',');
		write(v[i]);
	}
	putchar(']');
}
template<>
inline void write(const std::string&s){
	for(const auto&c:s)putchar(c);
}
template<typename T>
inline void writen(const T&f){
	write(f);putchar('\n');
}
template<typename T,typename...Args>
inline void writen(const T&f,const Args&...args){
	write(f);putchar(' ');writen(args...);
}
template<typename T>
inline void writeln(const T&f){
	write(f);putchar('\n');
}
template<typename T,typename...Args>
inline void writeln(const T&f,const Args&...args){
	write(f);writeln(args...);
}
template<typename T>
inline __int128 isread(T&t){t=0;bool f=0;char c=getchar();
	while(c!=EOF&&(c<'0'||c>'9'))f=c=='-',c=getchar();
	if(c==EOF)return -1;
	while(c!=EOF&&!(c<'0'||c>'9'))t*=10,t+=c-'0',c=getchar();
	ungetc(c,stdin);
	t=f?-t:t;
	return 1;
}
inline __int128 isread(std::string&s){
	s.clear();
	char c=getchar();
	while(c!=EOF&&(c==' '||c=='\t'||c=='\r'||c=='\n'))c=getchar();
	if(c==EOF)return -1;
	while(c!=EOF&&!(c==' '||c=='\t'||c=='\r'||c=='\n')){s+=c;c=getchar();}
	ungetc(c,stdin);
	return 1;
}
template<typename T1,typename T2>
inline __int128 isread(std::pair<T1,T2>&p){
	if(isread(p.first)==-1)return -1;
	return isread(p.second);
}
template<typename T>
inline __int128 isread(std::vector<T>&v){
	std::vector<T>t{};
	char c=getchar();
	while(c!=EOF&&(c==' '||c=='\t'||c=='\r'))c=getchar();
	bool bracket=0;
	if(c=='['){
		bracket=1;
		c=getchar();
	}
	while(true){
		if(bracket){
			while(c!=EOF&&c!='\n'&&c!=']'&&(c==' '||c=='\t'||c=='\r'||c==','||c=='('||c==')'))c=(char)getchar();
			if(c==EOF||c=='\n'||c==']')break;
		}else{
			while(c!=EOF&&c!='\n'&&(c==' '||c=='\t'||c=='\r'))c=getchar();
			if(c==EOF||c=='\n')break;
		}
		ungetc(c,stdin);
		T x{};
		if(isread(x)==-1)break;
		t.push_back(std::move(x));
		c=getchar();
	}
	if(bracket&&c==']'){
		char n=getchar();
		if(n!='\n')ungetc(n,stdin);
	}
	if(!bracket&&t.empty())return -1;
	v=std::move(t);
	return 1;
}
inline __int128 isreadn(){return 0;}
template<typename T,typename...Args>
inline __int128 isreadn(T&f,Args&...args){
	auto t1=isread(f),t2=isreadn(args...);
	return (t1==-1||t2==-1)?-1:t1+t2;
}
template<typename T=__int128,T N=26>
class AC{
	public:
	std::vector<std::array<T,N>>dt;
	std::vector<T>fail,endp,ans,ord;
	std::map<T,T>mp;
	AC():dt{},fail{}{dt.push_back({});}
	AC(const std::vector<std::array<T,N>>&_dt,const std::vector<T>_fail={}):dt{_dt},fail{_fail}{}
	void clear(){dt.assign(1,{});fail.clear();endp.clear();ans.clear();ord.clear();mp.clear();}
	AC&operator+=(const std::string&s){
		T u{};
		for(const auto&c:s){
			if(!dt[u][c-'a']){
				dt[u][c-'a']=dt.size();
				dt.push_back({});
			}
			u=dt[u][c-'a'];
		}
		mp[u]++;
		return *this;
	}
	AC&operator*=(const std::string&s){
		T u{};
		for(const auto&c:s){
			if(!dt[u][c-'a']){
				dt[u][c-'a']=dt.size();
				dt.push_back({});
			}
			u=dt[u][c-'a'];
		}
		endp.push_back(u);
		return *this;
	}
	void reload(){
		fail.assign(dt.size(),{});
		ans.assign(dt.size(),{});
		ord.clear();
		std::queue<T>q;
		for(T i=0;i<N;i++)
			if(dt[0][i])q.push(dt[0][i]);
		while(!q.empty()){
			auto u=q.front();q.pop();
			ord.push_back(u);
			for(T i=0;i<N;i++){
				if(dt[u][i]){
					fail[dt[u][i]]=dt[fail[u]][i];
					q.push(dt[u][i]);
				}
				else dt[u][i]=dt[fail[u]][i];
			}
		}
	}
	T operator()(const std::string&s){
		T u{},res{};
		for(const auto&c:s){
			u=dt[u][c-'a'];
			for(T p=u;p&&mp[p]!=-1;p=fail[p]){
				res+=mp[p];
				mp[p]=-1;
			}
		}
		return res;
	}
	T operator[](const std::string&s){
		T u{};
		for(const auto&c:s){
			u=dt[u][c-'a'];
			ans[u]++;
		}
		for(auto i=ord.size();i--;)
			ans[fail[ord[i]]]+=ans[ord[i]];
		T mx{};
		for(const auto&e:endp)
			if(ans[e]>mx)mx=ans[e];
		return mx;
	}
};
__int128 n;
std::string s;
AC ac;
int main(){
	while(1){
		isreadn(n);
		if(!n)return 0;
		ac.clear();
		std::vector<std::string>strs;
		while(n--){isreadn(s);strs.push_back(s);ac*=s;}
		ac.reload();
		isreadn(s);
		__int128 mx=ac[s];
		writeln(mx);
		for(__int128 i=0;i<ac.endp.size();i++)
			if(ac.ans[ac.endp[i]]==mx)
				writeln(strs[i]);
	}
	return 0;
}