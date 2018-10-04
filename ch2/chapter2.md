# 2. 信号与系统的时域分析
## 目录
### 第1章 信号与系统的基本概念 1
### 第2章 线性时不变系统的时域分析 1
### §2.0 引言 1
### §2.1 用单位冲激表示信号 1
### §2.2 离散时间LTI系统的时域分析 3
### §2.3 连续时间LTI系统的时域分析 5
### §2.4 单位冲激响应与LTI系统性质 6
### §2.5 LTI系统的微分、差分方程描述 9
### §2.6 LTI系统的零状态响应和零输入响应 12
### §2.7 由微分方程和差分方程藐视的LTI系统的方框图表示 16
## §2.0　引言
1)	上一章介绍了基本的信号，如冲激信号、复指数信号等，以及系统的不同性质
2)	我们研究的对象主要是线性时不变系统，LTI
-  A)	线性系统的特点是叠加性（可加性和其次性的组合）
-  B)	时不变就是随着输入信号的延迟，输出信号也延迟相应时间
3)	为能进行系统分析，将信号分解成基本信号的线性组合，如冲激信号的组合
4)	然后对于系统，我们就考查其冲激响应即可
5)	总之，本章目标：
-  A)	信号：用冲激信号的线性组合来表示
-  B)	系统：用冲激响应来表示
-  C)	然后学会一个信号通过一个系统之后会发生什么样的变化。
## §2.1　用单位冲激表示信号
1)	离散x[n]用冲激函数表示
- A)	考查一个离散信号x[n](下面补上图，奥本海默第一版P53 图3.1 (a))



- B)	上章:x[n]\delta[n]=x[0]\delta[n]=x[0]，和x[n]\delta[n-n_0]=x[n_0]\delta[n-n_0]=x[n_0],
### 上面信号的每一个分量都可以用该式来表达
- a) x[-1]\delta[n+1]=$$ f(n)= \begin{cases} x[-1], & \text {$n=-1$} \\ 0, & \text{$n\neq-1$} \end{cases} $$
- b) x[0]\delta[n]=$$ f(n)= \begin{cases} x[0], & \text {$n=0$} \\ 0, & \text{$n\neq0$} \end{cases} $$
- c)x[1]\delta[n-1]=$$ f(n)= \begin{cases} x[1], & \text {$n=1$} \\ 0, & \text{$n\neq1$} \end{cases} $$
- d)由此推广,可得: x[n]=\Lambda+x[-4]\delta[n+4]+x[-3]\delta[n+3]+x[-2]\delta[n+2]+x[-1]\delta[n+1]+x[0]\delata[n-1]+x[2]\delta[n-2]+x[3]\delta[n-3]+x[4]\delta[n-4]+\Lambda
- e)当n=0时,根据上式x[0]=x[0]\delta[n],其余各项均等于0.当n=1时，x[1]=x[1]\delta[n-1],说明上面表达式是正确的。
- C)所以有：$ x[n]= \sum_{k \to -\infty}^\infty \x[k]\delta[n-k] $这就意味着可以把任何一个序列表示成一串移位单位脉冲序列\delta[n-k]的线性组合，而在这个组合中，权因子就是x[k]
- a)	例，x[n]=u[n]，用上式表示，$ u[n]= \sum_{k \to -\infty}^\infty \u[k]\delta[n-k] $=$ x[n]= \sum_{k=0}^\infty \\delta[n-k] $
b)	将n－k用m表示，则 $ u[n]= \sum_{m=n}^\infty \\delta[m] $=$  \sum_{k \to -\infty}^\n \\delta[k] $，跟上一章的结果完全一样
2)	连续x(t)用冲激函数表示（两种方法，一种直观推导，一种用公式）
- A)	跟离散类似，连续时间函数也可用冲激函数来表示，形式也跟离散类似，不过是用积分表示：$$\int_{n \to -\infty}^{n \to +\infty} {x(t)\delta(t-\tau)} \,{\rm d}t$$ ，具体理解方式如下：
- B)	对于连续信号x（t），我们一般用\hat{x}(t) 来近似，图示
- a)	为了将\hat{x} 用封闭数学公式表示，我们引入以前的\delta_\Delta(t)=$ f(n)= \begin{cases} 1/\Delta, & \text {0<t<\delta} \\ 0, & \text{t^2为其他值} \end{cases} $$
面积为1的矩形，这样的话，就有\hat{x}(t)=\sum_{k \to -\infty}^\infty \x(k\Delta)\delta_\Delta(t-k\Delta)\Delta $
因为 ，\delta_\Delta(t)\Delta=1所以有上式
- b)	随着Δ \to 0，\hat{x}(t) 将愈来愈近似于x(t)，最后的极限是x(t)即有x(t)=