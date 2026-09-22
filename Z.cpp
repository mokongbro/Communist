#include<cstdio>
/**
 * @brief 输出一个整数（十进制，带符号）。
 *
 * @details 递归转换为十进制并逐位输出，顺序与人类阅读一致（高位在前）：
 * 1. \p t 为负时先 \c putchar('-')，再对相反数继续递归；
 * 2. \p t 小于 10 时直接输出这一位数字（\c t+'0'）；
 * 3. 否则先递归输出 \c t/10 的更高位，再输出最低位 \c t%10。
 * 全程只使用比较、除 10 与取模，对 \c __int128 等扩展整型同样适用；
 * 对 \c char 另有特化，按字符而非数值输出。
 *
 * @note 对类型的最小负值取相反数会溢出，属该写法已知的理论限制（正常数据不会触及）。
 * @tparam T 支持与 0 比较、除以 10 的整型（含 \c __int128 等扩展整型）。
 * @param t 待输出的值。
 */
template<typename T>
inline void write(const T&t){
	if(t<0)return putchar('-'),write(-t);
	if(t>9)write(t/10);
	putchar(t%10+'0');
}
/**
 * @brief \c char 特化：直接输出字符本身（而不是它的数值）。
 *
 * @details 泛型 \ref write 会把 \c char 视作整数按十进制输出——例如换行字符（数值 10）
 * 会被输出成两个字符 \c 10；本特化拦截该情形，直接 \c putchar 输出 \p c，
 * 使传入换行、空格等字符的 \c write 调用按字符原样输出。
 * \ref writeln 的变参重载正是依赖本特化来支持自定义分隔符的输出。
 * @param c 待输出的字符。
 */
template<>
inline void write(const char&c){
	putchar(c);
}
/**
 * @brief 输出一个值并换行。
 *
 * @details 等价于先 \ref write 输出 \p f，再由 \c putchar 补一个换行字符；
 * 变参重载递归展开的最终一步也会落到本重载，从而保证每行输出末尾恰好一个换行。
 * 需要在同一行输出多个值时，请改用变参重载或自行多次调用 \ref write。
 * @tparam T 任意可被 \ref write 输出的类型。
 * @param f 待输出的值。
 */
template<typename T>
inline void writeln(const T&f){
	write(f);putchar('\n');
}
/**
 * @brief 依次输出多个值（相互间无分隔符），最后换行。
 *
 * @details 递归展开：先 \ref write 输出第一个值 \p f，再对其余值递归调用本函数，
 * 直到只剩一个值落到单参重载，输出后补上换行。
 * 参数之间不会自动插入任何分隔符——需要分隔时把分隔字符本身作为一个参数传入，
 * 例如 \c writeln(ansz,'\\n',ansp) 会输出 ansz、换行、ansp、换行
 * （其中的换行字符依赖 \c write 的 \c char 特化按字符输出，而不是数值 10）。
 * @tparam T 第一个值的类型。
 * @tparam Args 其余值的类型。
 * @param f 第一个待输出的值。
 * @param args 其余待输出的值。
 */
template<typename T,typename...Args>
inline void writeln(const T&f,const Args&...args){
	write(f);writeln(args...);
}
#include<string>
#include<vector>
/**
 * @brief 读取一个整数。
 *
 * @details 全程 \c getchar 逐字符处理：
 * 1. 先清零 \p t，然后跳过所有非数字字符；跳过期间记录「最后一个非数字字符是否为负号」
 * 作为符号位（因此支持 \c "-123" 形式，行首空白与换行会被自然跳过）；
 * 2. 若跳过阶段直接到达 EOF，读入失败，返回 -1；
 * 3. 逐位累加 \c t=t*10+(c-'0')，直到非数字字符或 EOF；
 * 4. 把终止字符用 \c ungetc 放回输入流（为 EOF 时无副作用），最后按符号位取负。
 * 与字符串重载可混用：两者都会把终止字符放回，互不吞掉对方的输入。
 * @tparam T 目标整型（含 \c __int128）。
 * @param t 输出参数：读到的整数（按引用写回，先清零）。
 * @return 成功返回 1；在读到任何数字前遇到 EOF 返回 -1。
 */
template<typename T>
inline __int128 isread(T&t){t=0;bool f=0;char c=getchar();
	while(c!=EOF&&(c<'0'||c>'9'))f=c=='-',c=getchar();
	if(c==EOF)return -1;
	while(c!=EOF&&!(c<'0'||c>'9'))t*=10,t+=c-'0',c=getchar();
	ungetc(c,stdin);
	t=f?-t:t;
	return 1;
}
/**
 * @brief 读取一个以空白分隔的字符串。
 *
 * @details 即 \c cin>>s 的语义，用 \c getchar 实现：
 * 1. 先 \c clear 结果串，跳过所有空白字符（空格、制表符、回车、换行）；
 * 2. 若跳过阶段直接到达 EOF，读入失败，返回 -1；
 * 3. 用 \c s+=c 连续收集字符，直到遇到空白或 EOF 为止；
 * 4. 把终止字符用 \c ungetc 放回输入流，与后续读入（含 \ref isread 的其它重载）无缝衔接。
 * @note 空白字符本身不会被读入；读取只含小写字母的串时，效果等同按行读取。
 * @param s 输出参数：读到的字符串（先 \c clear）。
 * @return 成功返回 1；在读到任何字符前遇到 EOF 返回 -1。
 */
inline __int128 isread(std::string&s){
	s.clear();
	char c=getchar();
	while(c!=EOF&&(c==' '||c=='\t'||c=='\r'||c=='\n'))c=getchar();
	if(c==EOF)return -1;
	while(c!=EOF&&!(c==' '||c=='\t'||c=='\r'||c=='\n')){s+=c;c=getchar();}
	ungetc(c,stdin);
	return 1;
}
/**
 * @brief \ref isreadn 的变参递归终止重载。
 *
 * @details 供变参版本在「没有剩余参数」时调用，是递归展开的边界情形；
 * 返回 0 使成功计数之和从 0 开始累加。
 * @return 恒为 0，作为成功计数的初值。
 */
inline __int128 isreadn(){return 0;}
/**
 * @brief 对每个参数依次调用 \ref isread 读入多个值。
 *
 * @details 递归展开：先对第一个参数 \p f 调用 \ref isread，再对其余参数递归调用本函数；
 * 每个参数都会被尝试读取，最后只要有任一参数读入失败（返回 -1），整体即返回 -1，
 * 全部成功则返回成功次数之和（等于参数个数）。
 * 常用形式如 \c isreadn(n)、\c isreadn(a,b)、\c isreadn(n,m,x)，
 * 一次调用连续读入以空白分隔的多个值。
 * @tparam T 第一个参数的类型。
 * @tparam Args 其余参数的类型。
 * @param f 第一个待读取的对象。
 * @param args 其余待读取的对象。
 * @return 全部成功时返回成功次数之和（等于参数个数）；任一参数读到 EOF 返回 -1。
 */
template<typename T,typename...Args>
inline __int128 isreadn(T&f,Args&...args){
	auto t1=isread(f),t2=isreadn(args...);
	return (t1==-1||t2==-1)?-1:t1+t2;
}
/**
 * @brief 计算字符串 \p s 的 Z 函数。
 *
 * @details 定义 Z 函数 $z_i=\operatorname{LCP}(s,s[i..])$（$0\le i<|s|$，0 下标），
 * 即 \p s 与它从第 $i$ 位开始的后缀的最长公共前缀长度，并约定 $z_0=|s|$。
 * 实现按 Z-box（当前已匹配到的最右区间 $[\texttt{l},\texttt{r}]$）递推，
 * 每次转移均摊 $O(1)$，时间复杂度 $O(|s|)$，空间复杂度 $O(|s|)$。
 *
 * @tparam T 数组元素类型，默认 \c int。$z_i$ 的值不超过 $|s|$；
 * 大长度字符串请保持 4 字节宽度（换成 \c __int128 后内存与常数都会显著增大）。
 * @param s 待处理的字符串（按值传入）。
 * @return std::vector<T> 长度为 $|s|$ 的 Z 函数数组，满足 $\texttt{z[0]=|s|}$，
 * 且 $\texttt{z[i]}=\operatorname{LCP}(s,s[i..])$。
 *
 * @pre \p s 非空。
 * @note 参数按值传递；传入大串时可传右值（临时对象、拼接结果）省去一次拷贝。
 * @see exKMP()
 *
 * @par 示例
 * @code
 * auto z=Z_func("aaaabaa");//z 为 {7,3,2,1,0,2,1}
 * @endcode
 */
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
/**
 * @brief 扩展 KMP（exKMP）：求模式串 \p P 与文本串 \p S 每个后缀的最长公共前缀。
 *
 * @details 定义 $p_i=\operatorname{LCP}(P,S[i..])$（$0\le i<|S|$，0 下标），函数返回整个 $p$ 数组。
 * 实现把三个部分拼接为 $\texttt{P+'\backslash x7f'+S}$ 后一次性调用 Z_func()，
 * 再截取分隔符之后的 $|S|$ 个位置作为答案：分隔字符（ASCII 0x7F）与 $P,S$ 的字符均不同，
 * 匹配到它即停止，因此恒有 $p_i\le|P|$，结果与直接在 \p S 上跑 exKMP 等价。
 * 时间复杂度 $O(|S|+|P|)$，空间复杂度 $O(|S|+|P|)$。
 *
 * @tparam T 数组元素类型，默认 \c int（理由同 Z_func()）。
 * @param S 文本串（例如题目中的 $a$），按值传入。
 * @param P 模式串（例如题目中的 $b$），按值传入。
 * @return std::vector<T> 长度为 $|S|$ 的数组，$v_i=\operatorname{LCP}(P,S[i..])\le|P|$。
 *
 * @warning 要求 \p P 与 \p S 中不含 ASCII 0x7F 字符（小写字母串自然满足）。
 * @see Z_func()
 *
 * @par 示例
 * @code
 * auto p=exKMP("aaaabaa","aaaaa");//p 为 {4,3,2,1,0,2,1}
 * @endcode
 */
template<typename T=int>
std::vector<T>exKMP(std::string S,std::string P){
	T ns=(T)S.size(),nt=(T)P.size();
	std::vector<T>v(ns),z=Z_func(P+'\x7f'+S);
	for(T i=0;i<ns;i++)v[i]=z[nt+1+i];
	return v;
}
std::string a,b;///< 输入的两个字符串：a 为文本串，b 为模式串（各占一行）。
__int128 n,m,ansz,ansp;///< n=|a|，m=|b|；ansz、ansp 分别为 z 数组与 p 数组的权值（异或和，默认初值 0）。
/**
 * @brief 洛谷 P5410【模板】扩展 KMP / exKMP（Z 函数）解答的主程序。
 *
 * @details 读入两行字符串 \p a、\p b（$1\le|a|,|b|\le2\times10^7$，全为小写字母）后：
 * 1. 调用 Z_func(b) 求 $b$ 的 Z 函数数组 $z$；
 * 2. 调用 exKMP(a,b) 求 $p_i=\operatorname{LCP}(b,a[i..])$；
 * 3. 按定义计算权值：$z$ 的权值为 $\bigoplus_{i=1}^{m}i\cdot(z_i+1)$，
 * $p$ 的权值为 $\bigoplus_{i=1}^{n}i\cdot(p_i+1)$（$z_i,p_i$ 为 1 下标，
 * 对应代码中 0 下标的 $\texttt{(i+1)*(arr[i]+1)}$）；
 * 4. 依次输出 $z$ 的权值、$p$ 的权值，各占一行。
 * 总复杂度 $O(|a|+|b|)$；两个数组的元素类型为 \c int（4e7 规模下 \c __int128 数组会超出 500MB 内存限制）。
 *
 * @return 正常结束返回 0；读到文件结尾（\ref isreadn 返回 -1）时提前返回 0。
 * @see Z_func() exKMP() isreadn()
 */
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