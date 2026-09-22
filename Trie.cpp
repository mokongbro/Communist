#include<cstdio>
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
template<typename T>
inline void writeln(const T&f){
	write(f);putchar('\n');
}
template<typename T,typename...Args>
inline void writeln(const T&f,const Args&...args){
	write(f);writeln(args...);
}
#include<string>
#include<array>
#include<vector>
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
inline __int128 isreadn(){return 0;}
template<typename T,typename...Args>
inline __int128 isreadn(T&f,Args&...args){
	auto t1=isread(f),t2=isreadn(args...);
	return (t1==-1||t2==-1)?-1:t1+t2;
}
template<typename T=int,T N=62>
class Trie{
	public:
	std::vector<std::array<T,N>>dt;
	std::vector<T>cnt;
	Trie(){dt.push_back({}),cnt.push_back(0);}
	void clear(){dt.assign(1,{}),cnt.assign(1,0);}
	static int idx(const char&c){
		return c<='9'?c-'0'+52:c<='Z'?c-'A'+26:c-'a';
	}
	Trie&operator+=(const std::string&s){
		T u{};
		for(const auto&c:s){
			if(!dt[u][idx(c)]){
				dt[u][idx(c)]=(T)dt.size();
				dt.push_back({}),cnt.push_back(0);
			}
			u=dt[u][idx(c)],cnt[u]++;
		}
		return *this;
	}
	T operator()(const std::string&s)const{
		T u{};
		for(const auto&c:s){
			u=dt[u][idx(c)];
			if(!u)return 0;
		}
		return cnt[u];
	}
};
Trie tr;
__int128 T,n,q;
std::string s;
int main(){
	isread(T);
	while(T--){
		tr.clear();
		isreadn(n,q);
		while(n--)isread(s),tr+=s;
		while(q--)isread(s),writeln(tr(s));
	}
	return 0;
}