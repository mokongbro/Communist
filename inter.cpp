#include<map>
unsigned __int128 C(unsigned __int128 n,unsigned __int128 m){
	static std::map<std::pair<unsigned __int128,unsigned __int128>,unsigned __int128>M{};
	static unsigned __int128 len{};
	if(!M[{0,0}])M[{0,0}]=1;
	if(M[{n,m}])return M[{n,m}];
	while(len<n){
		len++;
		M[{len,0}]=1;
		for(unsigned __int128 i=1;i<=len;i++)M[{len,i}]=M[{len-1,i-1}]+M[{len-1,i}];
	}
	return M[{n,m}];
}