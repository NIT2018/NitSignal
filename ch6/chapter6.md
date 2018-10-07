
目       录
第6章 拉普拉斯变换	1
§6.0 引言	1
§6.0.1 傅立叶变换概述	1
§6.0.2 本章目标	2
§6.1 拉普拉斯变换	2
§6.2 拉普拉斯变换的收敛域(ROC: Region of Convergence)	3
§6.3 收敛域的相关性质	3
§6.4 拉普拉斯反变换	7
§6.4.1 拉普拉斯反变换公式	7
§6.4.2 常用变换对	7
§6.4.3 几个例子	7
§6.5 拉普拉斯变换的性质	8
§6.5.1 拉氏变换的线性性质	8
§6.5.2 时域平移	8
§6.5.3 s域平移	8
§6.5.4 时间尺度变换	9
§6.5.5 卷积性质	9
§6.5.6 时域微分	9
§6.5.7 s域微分	10
§6.5.8 时域积分	10
§6.5.9 初值与终值定理	10
§6.6 利用拉氏变换分析和表征LTI系统	11
§6.6.1 系统函数	11
§6.6.2 由线性常系数微分方程表征的系统	11
§6.6.3 LTI系统互联后的系统函数	12
§6.6.4 用拉氏变换进行系统分析	13
§6.7 单边拉氏变换	13


# 第6章　拉普拉斯变换
## §6.0　引言
### §6.0.1　傅立叶变换概述
1)在傅立叶变换中，我们已经知${\it e}^{st} \to {\mit H}(s){\it e}^{st}$，其中${\it e}^{st}$是系统特征函数，H(s)则是特征值；并且s=jw；

2)拉普拉斯将傅立叶的结论进行了扩展，$s=\sigma+j\omega $，这样就成为了拉普拉斯变换

3)因此，两者在很大程度上有相似性
### §6.0.2　本章目标
1)掌握拉普拉斯变换，性质，定义域

2)掌握拉普拉斯变换在LTI系统分析中的应用
### §6.1　拉普拉斯变换
1)根据特征值和特征函数的关系，有${\it e}^{st} \to {\mit H}(s){\it e}^{st}$，其中$\int_{+\infty}^{-\infty} h(t)e^{st} \{\rm d}t$，就是h(t)的拉普拉斯变换
2)所以x(t)的拉普拉斯定义：$X(s)=\Delta\int_{+\infty}^{-\infty}{\it x}(t){\it e}^{-st}\{\rm d}t$，记为$x(t)\overleftrightarrow{f} X(s)$
3)当s=jw时，就变成了傅立叶变换，所以傅立叶变换是拉普拉斯的一种特例

4)例1：求信号$x(t)={\it e}^{-at}u(t),{\it a}>0$的拉普拉斯变换

a)我们已知该信号的傅立叶变换为$X(j\omega)=\frac{1}{a+j\omega}$
i:  $$ X(s)=\int_{+\infty}^{-\infty} {x(t) e^{-st}dt = \int_{+\infty}^{0}e^{-at} e^{(-\sigma - j\omega) t}dt  =  \frac{{1}{-a-\sigma-j\omega}} \left. e^{(-a-\sigma-j\omega)t} \right| _{+\infty}^{x=0} = \frac{1}{a + \sigma - j\omega} = \frac{1}{s+a} $$
，
其中必须满足条件$\sigma+a>0，\sigma>-a,Re{\it s}>-a$,即

b)傅立叶变换是Re{s}=0的情况
5)例2：求信号$x(t)=-{\it e}^{-at}u(-t),a>0$的拉普拉斯变换

a)解：$$
，
其中必须满足条件$\sigma+a<0,\sigma<-a,Re{\it s}<-a$

6)从以上两个例子我们看出，同样的拉普拉斯，有可能是不同的原信号，这里有一个非常重要的收敛域的概念

7)同时，如果虚轴Re{s}=0包含在收敛域内的话，则对应信号的傅立叶变换存在，否则不存在
## §6.2　拉普拉斯变换的收敛域(ROC: Region of Convergence)
1)$X(s)=\Delta\int_{+\infty}^{-\infty}{\it x}(t){\it e}^{-st}\{\rm d}t$

2)能使拉普拉斯变换收敛的s值的范围称为拉普拉斯变换的收敛域，一般简称为ROC

3)一般来说，要使其收敛主要考察s的实部即可
## §6.3　收敛域的相关性质
1)性质1：X(s)的ROC在s平面内由平行于jω轴的带状区域所组成

a)因为X(s)的ROC是由这样一些s=σ+jω所组成，在这些范围里面，要使$x(t){\it e}^{-\sigma t}$的傅立叶变换收敛，因此ROC只跟s的实部有关，所以ROC是平行于jω轴的带状区域

2)性质2：对有理拉氏变换来说，ROC内部不包括任何极点

a)有理拉氏变换概念。，$X(s)=\frac{N(s)}{D(s)}$,N(s)和D(s)均是s的多项式，X(s)具有这种形式的时候，叫有理拉氏变换。N(s)=0的根叫做零点，D(s)=0的根叫做极点。

b)所以我们可以理解ROC内部不能包含极点，否则的话，在该点处的X(s)值变成无限大，显然不收敛。

3)性质3：如果x(t)是有限持续期，并至少存在一个s，使其拉氏变换收敛，那么ROC就是整个s平面

a)x(t)是有限持续期，也就是在[T1，T2]，期间，x(t)取有限值，而在这个范围之外，x(t)=0，如图所示(奥本海默第一版P443图9.4)：

b)现在假设对于某个s=σ0+jω来说，拉氏变换收敛，也就是$x(t){\it e}^{-\sigma_0 t}$的傅立叶变换存在，也即$x(t){\it e}^{-\sigma_0 t}$绝对收敛。所以有$\frac{T_2}{T_1}{|x(t)|}{\it e}^{-\sigma_0 t}\{\rm d}t<\infty$

c)接下来，我们证明对于$\sigma_1>\sigma_0$，或者$\sigma_1<\sigma_0$，均能使$x(t){\it e}^{-\sigma_1 t}$绝对收敛，这样就可以证明X(s)的收敛域是整个s平面。

d)先假设$\sigma_1>sigma_0$，则$$\int_{T_2}^{T_1}{|x(t)|}{\it e}^{-\sigma_1 t}\{\rm d}t=\int_{T_2}^{T_1}{|x(t)|}{\it e}^{-\sigma_0 t}\{\it e}^{-(\sigma_1-\sigma_0)t}\{\rm d}t<{\it e}^{-(\sigma_1-\sigma_0)T_1}\int_{T_2}^{T_1}{|x(t)|}{\it e}^{-\sigma_0 t}\{\rm d}t<\infty$$
，所以$x(t){\it e}^{-\sigma_1 t}$绝对收敛，因此对于$\sigma_1>\sigma_0$，如果在s=σ0+jω收敛，那么在s=σ1+jω处也收敛

e)同理，当$\sigma_1<\sigma_0$时，也可以证明$x(t){\it e}^{-\sigma_1 t}$收敛。$$\int_{T_2}^{T_1}{|x(t)|}{\it e}^{-\sigma_1 t}\{\rm d}t=\int_{T_2}^{T_1}{|x(t)|}{\it e}^{-\sigma_0 t}\{\it e}^{-(\sigma_1-\sigma_0)t}\{\rm d}t<{\it e}^{-(\sigma_1-\sigma_0)T_1}\int_{T_2}^{T_1}{|x(t)|}{\it e}^{-\sigma_0 t}\{\rm d}t<\infty$$

f)所以原命题成立。

g)对于该命题，我们需要强调的是两个条件，1：x(t)是有限持续期；2：至少存在一个s值，使其拉氏变换收敛

4)性质4：如果x(t)是右边信号，而且如果Re{s}=σ0这条线位于ROC内，那么Re{s}>σ0的全部s都一定在ROC内。

a)若在某个时刻T1以前，x(t)=0，那么该信号就称为右边信号，如图所示（奥本海默第一版P444图9.6）

b)同样，我们假设对于某个s=σ0+jω来说，拉氏变换收敛，也就是$x(t){\it e}^{-\sigma_0 t}$的傅立叶变换存在，也即$x(t){\it e}^{-\sigma_0 t}$绝对收敛。所以有$\int_{\infty}^{T1}\{|x(t)|}{\it e}^{-\sigma_0 t}\{\rm d}t<\infty$

c)对于$\sigma_1>\sigma_0$，我们也要证明$x(t){\it e}^{-\sigma_1 t}$绝对收敛。$$\int_{+\infty}^{T_1}{|x(t)|}\{\it e}^{-\sigma_1 t}\{\rm d}t=\int_{+\infty}^{T_1}{|x(t)|}{\it e}^{-\sigma_0 t}\{\it e}^{-(\sigma_1-\sigma_0)t}\{\rm d}t<{\it e}^{-(\sigma_1-\sigma_0)T_1}\int_{+\infty}^{T_1}{|x(t)|}{\it e}^{-\sigma_0 t}\{\rm d}t<\infty$$

d)所以原命题成立。

5)性质5：如果x(t)是左边信号，而且如果Re{s}=σ0这条线位于ROC内，那么Re{s}<σ0的全部s都一定在ROC内。

a)若在某个时刻T1以后，x(t)=0，那么该信号就称为左边信号，如图所示（奥本海默第一版P445图9.8）

b)完全类同性质4，此处从略。

6)性质6：如果x(t)是双边信号，而且如果Re{s}=σ0这条线位于ROC内，那么ROC就一定是由s平面的一条带状区域所组成，直线Re{s}=σ0位于带中。

a)一个双边信号就是对t>0和t<0都具有无限范围的信号，如图所示（奥本海默第一版P445图9.9(a)）。

b)这种信号的ROC如何求呢？

c)我们取个时间点T0，把信号分成t>T0和t<T0两个叠加，一个是右边信号，另一个左边信号。分别应用以上性质4和性质5。

d)若s=σL+jω是左边信号的一个收敛范围，那么σ<σL所有范围均是收敛域范围；同理若s=σR+jω是右边信号的一个收敛范围，那么σ>σR所有范围均是收敛域范围。若σL>σR，则有共同区域（是一个带状区域），否则没有(当然，这里的σL和σR也代表了收敛域的边界值)，如图所示

e)所以原命题是成立

7)例1：求信号$$ f(n)= \begin{cases} {\it e}^{-at}, & \text {0<t<T} \\ 0, & \text{其他t值} \end{cases} $$对应拉普拉斯变换X(s)的零极点

a)解：$$X(s)=\int_{+\infty}^{-\infty}
x(t){\it e}^{-st}\{\rm d}t=\int_{T}^{0}\{\it e}^{-at}\{\it e}^{-st}\{\rm d}t=frac{1}{s+a}\(1-{\it e}^{-(s+a)T})$$

b)再确定ROC。根据性质2，有限持续的信号，其ROC会覆盖全部s平面。不过从上面的X(s)表达式中，我们似乎看到s=-a是一个极点，这个跟ROC中不包含极点相冲突。其实不然。

c)当s=-a时，分母为0，分子也为零。当s=-a时，X(s)的值可以通过求s\rightarrow-a时的极限来得到。$$ \lim_{s \to -a}X(s)=\lim_{s\to -a}\frac[{\frac{d}{ds}\(1-{\it e}^{-(s+a)T})}{\frac{d}{ds}\(s+a)}]=\lim_{s\to -a}T{\it e}^{-(s+a)T}=T $$，所以X(-a)=T，所以s=-a既不是零点也不是极点。

d)接下来求零点。零点满足条件为：${\it e}^{-(s+a)T}=1={\it e}^{-j2k\pi}$，也即，$(s+a)T=j2k\pi$,$s=-a+j\frac{2k\pi}{T},k=0,\pm 1,\pm 2\ldots$

e)根据上面分析s=-a不是零点，所以，零点为：$s=-a+j\frac{2k\pi}{T},k=0,\pm 1,\pm 2\ldots$
，如下图所示。

8)例2：求信号$x(t)={\it e}^{-b|t|}$的拉普拉斯变换

a)解：将x(t)分解成左边信号和右边信号之和。$$x(t)={\it e}^{-b|t|}={\it e}^{bt}u(t)+{\it e}^{bt}u(-t)$$

b)由前面介绍，有${\it e}^{-bt}u(t)\leftrightarrow \frac{1}{s+b}$,Re{s}<-b

c)同理，${\it e}^{bt}u(t)\leftrightarrow \frac{-1}{s-b}$,Re{s}<b


d)所以有：${\it e}^{-b|t|}\leftrightarrow\frac{1}{s+b}\-\frac{1}{s-b}=\frac{2b}{s^2 -b^2}$,-b<Re{s}<b

e)因此，当b>0时，x(t)的拉普拉斯变换存在，公式和ROC如上所示。当b<0时，拉普拉斯变换不存在。

9)例3：已知信号x(t)的拉普拉斯变换为$X(s)=\frac{1}{(s+1)(s+2)}$，求原信号x(t)

a)解：信号的拉普拉斯变换，除了表达式之外，还要跟ROC结合起来。$X(s)=\frac{1}{s+1}\-\frac{1}{s+2}$，再考查X(s)的ROC。

b)由已知，信号有两个极点s=-1，s=-2，又根据性质，极点不能在ROC区域内，所以，有下列三种不同情况的ROC

<<<<<<< HEAD
c)情况(a)，右边信号，可得$x(t)={\it e}^{-t}u(t)-{\it e}^{-2t}u(t)$

d)情况(b)，左边信号，可得$x(t)={\it -e}^{-t}u(-t)+{\it e}^{-2t}u(-t)$

e)情况(c)，双边信号，可得$x(t)={\it -e}^{-t}u(-t)-{\it e}^{-2t}u(t)$


## §6.4　拉普拉斯反变换
### §6.4.1　拉普拉斯反变换公式
1)由前面介绍，我们知道$s=\sigma+j\omega$，要求拉普拉斯反变换的话，可以从傅立叶反变换的角度来考虑；

2)现在有$$X(s)=X(\sigma+j\omega)=F{X(t){\it e}^{-\sigma t}}=\int_{+\infty}^{-\infty}x(t){\it e}^{-\sigma t}{\it e}^{-j\omega t}\{\rm d}t$$

3)所以$$x(t){\it e}^{-\sigma t}=F^{-1}{X(\sigma +j\omega)}=\frac{1}{2\pi}\int_{+\infty}^{-\infty}X(\sigma+j\omega){\it e}^{j\omega t}\{\rm d}\omega$$

4)也即$$x(t)=\int_{+\infty}^{-\infty}X(\sigma +)j\omega){\it e}^{(\sigma +)j\omega)t}\{\rm d}\omega=\frac{1}{2\pi j}\int_{\sigma+j\infty}^{\sigma-j\infty}X(s){\it e}^{st}\{\rm d}s$$

5)这个积分的求值，是要求利用复平面的围线积分

6)一般求拉普拉斯反变换的话，可以利用常用的变换对
### §6.4.2　常用变换对
1)课本P235
### §6.4.3　几个例子
待完成课后习题后补充。
## §6.5　拉普拉斯变换的性质
### §6.5.1　拉氏变换的线性性质
1)若$x_1(t)\leftrightarrow X_1(s)$，ROC为R1；$x_2(t)\leftrightarrow X_2(s)$，ROC为R2。
则有$ax_1(t)+bx_2(t)\leftrightarrow aX_1(s)+bX_2(s)$，ROC为$R_1\bigcap R_2$
2)利用拉氏变换的公式$X(s)=\int_{+\infty}^{-\infty}x(t){\it e}^{-st}\{\rm d}t$即可证明。
3)如果有零极点抵消情况发生，则有可能最终的相交后ROC会比原来的ROC要大。
4)例：$$X_1(s)=\frac{1}{s+1}；Re{s}>-1；X_2(s)=\frac{1}{(s+1)(s+2)},Re{s}>-1；x(t)=x_1(t)-x_2(t)$$求X(s)
5)解：由拉氏变换的线性性质，得
.$$X(s)=X_1(s)-X_2(s)=\frac{1}{s+1}-\frac{1}{(s+1)(s+2)}=\frac{1}{s+2},Re{s}>-2$$
6)原来两个ROC的交集为Re{s}>-1，现在因为零极点抵消，所以ROC扩展到Re{s}>-2
### §6.5.2　时域平移
1)若$x(t)\leftrightarrow X(s)$，ROC=R；则有$x(t-t_0)\leftrightarrow {\it e}^{-st_0}X(s)$，ROC=R
### §6.5.3　s域平移
1)若$x(t)\leftrightarrow X(s)$，ROC=R；则有${\it e}^{s_0 t}x(t)\leftrightarrow X(s-s_0)$，ROC：R1＝R+Re{s0}
2)ROC的变化如下图所示。

### §6.5.4　时间尺度变换
1)若$x(t)\leftrightarrow X(s)$，ROC=R；则有$x(at)\leftrightarrow \frac{1}{|a|}X(\frac{s}{a})$，ROC：R1＝aR
2)ROC的变化如下图所示。

### §6.5.5　卷积性质
1)若$x_1(t)\leftrightarrow X_1(s)$，ROC为R1；$x_2(t)\leftrightarrow X_2(s)$，ROC为R2。则有$x_1(t)\ast x_2(t)\leftrightarrow X_1(s)X_2(s)$，ROC为$R_1\bigcap R_2$
2)如果有零极点抵消情况发生，则有可能最终的相交后ROC会比原来的ROC要大。
3)例：$X_1(s)=\frac{s+1}{s+2},Re{s}>-2；X_2(s)=\frac{s+2}{s+1},Re{s}>-1$；求$x(t)=x_1(t)\ast x_2(t)$的频谱；
4)解：由卷积性质得：$X(s)=X_1(t)X_2(t)=1$；其ROC为整个s平面；
### §6.5.6　时域微分
1)若$x(t)\leftrightarrow X(s)$，ROC=R；则有$\frac{dx(t)}{dt}\leftrightarrow sX(s)$，ROC：包括R
2)证明：$x(t)=\frac{1}{2\pi j}\int_{\sigma+j\omega}^{\sigma-j\omega}X(s){\it e}^{st}\{\rm d}s$，两边对t求导，即可得：
3)$\frac{dx(t)}{dt}=\frac{1}{2\pi j}\int_{\sigma+j\omega}^{\sigma-j\omega}X(s){\it e}^{st}\{\rm d}s$，也即$\frac{dx(t)}{dt}\leftrightarrow sX(s)$
4)新的ROC可能区域会更大
### §6.5.7　s域微分
1)若$x(t)\leftrightarrow X(s)$，ROC=R；则有$-tx(t)\leftrightarrow \frac{dX(s)}{ds}$，ROC：包括R
2)证明：$X(s)=\int_{+\infty}^{-\infty}x(t){\it e}^{-st}\{\rm d}t$，两边对s求导，得：$\frac{dX(s)}{ds}=\int_{+\infty}{-\infty}(-t)x(t){\it e}^{-st}\{\rm d}t$，也即$-tx(t)\leftrightarrow \frac{dX(s)}{ds}$
3)例1：求下列信号的拉氏变换。$x(t)=t{\it e}^{-at}u(t)$
4)解，由前面可知：${\it e}^{-at}u(t)\leftrightarrow \frac{1}{s+a}$,Re{s}>-a$；再由上面性质，可得：$$t{\it e}^{-at}u(t)\leftrightarrow \frac{-d}{ds}(\frac{1}{s+a})=\frac{1}{(s+a)^2}$$,Re{s}>-a

5)再用一次上面的性质，可得：$\frac{t^2}{2}{\it e}^{-at}u(t)\leftrightarrow\frac{1}{(s+a)^3}$,Re{s}>-a
6)更一般的有：$\frac{t^{(n-1)}{(n-1)!}{\it e}^{-at}u(t)\leftrightarrow\frac{1}{(s+a)^n}$,Re{s}>-a
### §6.5.8　时域积分
1)若$x(t)\leftrightarrow X(s)$，ROC=R；则有$\int_{t}^{-\infty}x(\tau)d\tau\leftrightarrow\frac{1}{s}X(s)$，ROC：$R\bigcap{Re{s}>0}$
2)证明：$u(t)\leftrightarrow\frac{1}{s}$，$x(t)\ast u(t)=\int_{t}{-\infty}x(\tau)d\tau$，所以利用卷积性质，可得：$\int_{t}{-\infty}x(\tau)d\tau\leftrightarrow\frac{1}{s}X(s)$

### §6.5.9　初值与终值定理
1)若t<0时，有x(t)=0。而且在t=0时刻，x(t)不包含任何冲激或者高阶奇异函数，那么就可以直接从拉氏变换中计算初值x(0+)[也就是x(t)当t从正值方向趋向于0时的值]和终值，即$t\rightarrow\infty$的x(t)值。
2)初值定理：$x(0^+)=\lim_{s\to \infty}sX(s)$；
3)终值定理：$\lim{t\to \infty}x(t)=\lim_{s\to 0}sX(s)$
## §6.6　利用拉氏变换分析和表征LTI系统
### §6.6.1　系统函数
1)根据傅立叶变换，我们知道系统输出的频谱是输入的频谱和系统频响的乘积。也即：
2)相对应，我们有，当s＝jw时，就是傅立叶变换的公式。此处，H(s)就系统函数（或转移函数），（而H(jw)叫系统频率响应）
3)为了能够唯一得到h(t)，则还需要说明H(s)的ROC，否则，只知道H(s)而不知道ROC的话，则无法确知系统
4)一般为了能够确定H(s)对应的ROC，我们给系统加上另外性质，如因果系统，或者稳定系统，或者因果稳定系统，这样就可根据这些信息来得到h(t)。
5)如一个因果稳定系统的系统函数为，，求系统冲激响应。
a)解，由已知，系统是因果稳定的，我们可知：因果，当t<0时，h(t)＝0，h(t)是右边信号；稳定：H(s)包含s平面虚轴；
b)再利用拉普拉斯变换的特性，可知该H(s)只有两个可能的ROC：Re{s}>-1或Re{s}<-1。
c)显然，该系统函数H(s)对应的ROC应该是Re{s}>-1，所以得到对应的h(t)为：
§6.6.2　由线性常系数微分方程表征的系统
1)基本思路跟傅立叶分析类似
2)例：求因果系统的系统函数H(s)。
a)解：两边进行拉普拉斯变换，可得：，因为是因果系统，所以h(t)是右边信号，所以对应H(s)的ROC应该是最右边极点的右边区域，也即ROC：Re{s}>-3
3)更一般情况，一个系统的方程为，其对应的系统函数为：。
4)零点就是方程的解。
5)极点就是方程的解。
§6.6.3　LTI系统互联后的系统函数
1)系统并联；如下图所示。

有，则。
2)系统级联；

有，则。
3)反馈联接；

a)有，，即，即，得
b)另法：，
，
所以，
§6.6.4　用拉氏变换进行系统分析
此处用几个课后作业题来说明。
§6.7　单边拉氏变换
1)一般来说，信号x(t)是单边信号，且在t<0时，有x(t)=0，因此，可以用到拉氏变换来处理；
2)单边拉氏变换公式：
3)此处的积分下限0，可以是0-，也可以是0+。一般为了考虑x(t)在原点处有冲激函数及其各阶导数，我们取积分范围为0-，这个在分析具有非零初始条件的LTI系统中有带来较大方便。所以实际上是，记做
4)例1：求的单边和双边拉氏变换；
a)解：双边拉氏变换：
，ROC：Re{s}>-a
b)单边拉氏变换：，ROC：Re{s}>-a
c)两者一样。
5)例2：求的单边和双边拉氏变换；具体解法
a)解：双边拉氏变换：
，ROC：Re{s}>-a
b)解单边拉氏变换：，ROC：Re{s}>-a
6)总结：对于t<0时，x(t)=0的信号，单边和双边拉氏变换是一样的，否则，就不一样。
7)时域微分性质：若，则
a)证明：

b)同理可得：
c)证明：

8)时域卷积性质：如x1(t)和x2(t)都是单边信号，即当t<0时，x1(t)=x2(t)=0。有，这里一定要注意公式的适用范围是x1(t)和x2(t)都是单边信号。
9)在实际分析中，对于都是单边信号的x(t)，都可以用单边拉氏变换来分析。