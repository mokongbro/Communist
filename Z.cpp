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
template<typename T=int>
std::vector<T>Z_func(std::string s){
	T n=(T)s.size();
	std::vector<T>z(n);z[0]=n;
	for(T i=1,l=0,r=0;i<n;i++){
		if(i<=r)z[i]=z[i-l]<r-i+1?z[i-l]:r-i+1;
		while(i+z[i]<n&&s[z[i]]==s[i+z[i]])z[i]++;
		if(i+z[i]-1>r)l=i,r=i+z[i]-1;
	}
	return z;
}
template<typename T=int>
std::vector<T>exKMP(std::string S,std::string P){
	T ns=(T)S.size(),nt=(T)P.size();
	std::vector<T>v(ns),z=Z_func(P+'\x7f'+S);
	for(T i=0;i<ns;i++)v[i]=z[nt+1+i];
	return v;
}
std::string a,b;
__int128 n,m,ansz,ansp;
int main(){
	if(isreadn(a,b)==-1)return 0;
	n=a.size(),m=b.size();
	auto p=exKMP(a,b);
	auto z=Z_func(b);
	for(__int128 i=0;i<m;i++)ansz^=(i+1)*(z[i]+1);
	for(__int128 i=0;i<n;i++)ansp^=(i+1)*(p[i]+1);
	writeln(ansz,'\n',ansp);
	return 0;
}