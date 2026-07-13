#include<map>
#include<vector>
#include<cmath>
unsigned __int128 C(unsigned __int128 n,unsigned __int128 m){
	static std::map<std::pair<unsigned __int128,unsigned __int128>,unsigned __int128>M{};
	static unsigned __int128 len{};
	if(m>n)return 0;
	if(n==0||m==0||m==n)return 1;
	if(M[{n,m}])return M[{n,m}];
	if(!M[{0,0}])M[{0,0}]=1;
	while(len<n){
		len++;
		M[{len,0}]=1;
		M[{len,len}]=1;
		for(unsigned __int128 i=1;i<len;i++)M[{len,i}]=M[{len-1,i-1}]+M[{len-1,i}];
	}
	return M[{n,m}];
}
extern "C"{
	extern long double Deri(long double f(long double),unsigned __int128 n,long double x=0){
		long double dx=(1+fabsl(x))*powl(1e-19l,1.0l/(n+1));
		static std::vector<long double>Pow;
		if(Pow.empty()||Pow[1]!=dx){
			Pow.clear();
			Pow.push_back(1);
			Pow.push_back(dx);
		}
		while(Pow.size()<=static_cast<size_t>(n))Pow.push_back(dx*Pow.back());
		long double ans{};
		for(unsigned __int128 i=0;i<=n;i++){
			unsigned __int128 coeff=C(n,i);
			long double term=f(x+static_cast<long double>(i)*dx);
			if((n-i)&1)ans-=static_cast<long double>(coeff)*term;
			else ans+=static_cast<long double>(coeff)*term;
		}
		return ans/Pow[static_cast<size_t>(n)];
	}
	extern long double inter(long double f(long double),long double a,long double b){
		const unsigned long long N=1000000ull*static_cast<unsigned long long>(fabsl(b-a));
		unsigned long long N_even=N;
		if(N_even%2)N_even++;
		long double h=(b-a)/N_even;
		long double sum=f(a)+f(b);
		for(unsigned long long i=1;i<N_even;i++){
			long double x=a+i*h;
			if(i&1)sum+=4.0l*f(x);
			else sum+=2.0l*f(x);
		}
		return sum*h/3.0l;
	}
}