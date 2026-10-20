#include<cstdio>
#include<cstdlib>
int main(int argc,const char*argv[]){
	if(argc!=5)
		for(int i=0;i<argc;i++)printf("%s\n",argv[i]);
	else{
		freopen(argv[1],"r",stdin);
		freopen(argv[2],"w",stdout);
		int inc=std::atoi(argv[3]),trc=std::atoi(argv[4]);
		int c=getchar();
		while(c!=EOF){
			putchar(c==inc?trc:c);
			c=getchar();
		}
	}
	return 0;
}