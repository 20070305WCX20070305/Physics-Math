# Wigner 因果界

> 主题：量子散射中的因果性约束。本文从薛定谔方程与量子力学基本原理出发，详细推导 Wigner 因果界，并解释其物理意义。
>
> 参考：E. P. Wigner, Phys. Rev. **98**, 145 (1955)；H.-W. Hammer and D. Lee, Phys. Lett. B **681**, 500 (2009)。

---

## 目录

1. 引言与大意
2. 散射理论与相移（从薛定谔方程出发）
3. 时间延迟：Eisenbud–Wigner 公式
4. 因果性原理
5. 倒数对数导数的单调性（Wronskian 严格推导）
6. Wigner 因果界的推导
7. 低能极限与有效程的 Wigner 界
8. 用 Bethe 公式交叉验证
9. 推广到任意维度 $d$ 与角动量 $L$
10. 物理意义与应用

---

## 1. 引言与大意

在量子力学中，**因果性**要求：入射波尚未到达散射体之前，不能有散射波传播出来。这一看似朴素的约束，对散射振幅、相移乃至低能有效程都施加了严格的界。

Eisenbud 和 Wigner 首先把"散射时间延迟"与相移对能量的导数联系起来：

$$
\Delta t=2\hbar\frac{d\delta}{dE}.
$$

因果性要求这个延迟不能是任意大的负值（时间超前）：入射波必须先进入相互作用区域，散射波才能出来。由于相互作用有有限程 $R$，最大时间超前被"穿越相互作用区所需时间"限制。由此得到相移导数的下界（**Wigner 因果界**）：

$$
\frac{d\delta_0}{dk}\ge -R+\frac{1}{2k}\sin\!\big[2(kR+\delta_0)\big],
$$

其在零能极限下化为对有效程 $r_0$ 的界：

$$
\boxed{\;r_0\le 2R-\frac{2R^2}{a}+\frac{2R^3}{3a^2}\;}
$$

其中 $a$ 为散射长度、$R$ 为势的程。本文逐步推导上述全部结果。

---

## 2. 散射理论与相移（从薛定谔方程出发）

### 2.1 径向薛定谔方程

考虑质量为 $m$ 的粒子在球对称势 $V(r)$ 中散射，$V(r)$ 为有限程势（当 $r>R$ 时 $V(r)=0$）。定态薛定谔方程为

$$
-\frac{\hbar^2}{2m}\nabla^2\psi+V(r)\psi=E\psi,\qquad E=\frac{\hbar^2k^2}{2m}.
$$

分离变量 $\psi(\mathbf r)=\dfrac{u_l(r)}{r}Y_{lm}(\theta,\varphi)$。利用 $\nabla^2\dfrac{u}{r}Y_{lm}=\left[\dfrac{1}{r}\dfrac{d^2u}{dr^2}-\dfrac{l(l+1)}{r^3}u\right]Y_{lm}$（因 $\nabla^2$ 径向部分作用在 $u/r$ 上给出 $\dfrac{1}{r}\dfrac{d^2}{dr^2}u$），得到径向方程

$$
\boxed{\;\left[-\frac{\hbar^2}{2m}\frac{d^2}{dr^2}+\frac{\hbar^2l(l+1)}{2mr^2}+V(r)\right]u_l(r)=E\,u_l(r)\;}
$$

边界条件 $u_l(0)=0$（波函数在原点正则）。

### 2.2 分波展开

入射平面波可展开为分波：

$$
e^{ikz}=\sum_{l=0}^{\infty}(2l+1)i^l j_l(kr)P_l(\cos\theta),
$$

其中 $j_l$ 为球贝塞尔函数，$P_l$ 为勒让德多项式。自由径向函数 $r j_l(kr)$ 在大 $r$ 处渐近为

$$
r j_l(kr)\sim\sin\!\left(kr-\frac{l\pi}{2}\right).
$$

### 2.3 相移

相互作用使第 $l$ 分波在远处多出一个相位 $\delta_l$：

$$
\boxed{\;u_l(r)\sim\sin\!\left(kr-\frac{l\pi}{2}+\delta_l\right)\qquad(r\to\infty)\;}
$$

$\delta_l=\delta_l(k)$ 称为**散射相移**。对低能散射，$l=0$（$s$ 波）占主导，本文以下主要讨论 $s$ 波，角标 $0$ 保留以便与一般情形对照。

### 2.4 散射振幅与 $S$ 矩阵

$s$ 波散射振幅为

$$
f_0(k)=\frac{e^{2i\delta_0}-1}{2ik}=\frac{1}{k\cot\delta_0-ik}.
$$

$S$ 矩阵元为

$$
S_0(k)=e^{2i\delta_0(k)}.
$$

注意：$S$ 矩阵携带的相位是 $2\delta_0$（出射波相对入射波），这一点在时间延迟的推导中至关重要。

---

## 3. 时间延迟：Eisenbud–Wigner 公式

### 3.1 波包的构造

构造一个由散射分波组成的波包。散射出射波（大 $r$）可写为

$$
\psi_{\rm sc}(r,t)=\frac1r\int dk\,g(k)\,S_0(k)\,e^{i(kr-\omega t)},
$$

其中 $g(k)$ 为波包权重（在 $k_0$ 附近尖锐），$S_0(k)=e^{2i\delta_0(k)}$，$\omega(k)=\hbar k^2/2m$。

### 3.2 稳相法求波包中心

相位为

$$
\Phi(k)=kr+2\delta_0(k)-\omega(k)t.
$$

波包在 $r$ 处出现（波包中心）的条件是相位平稳：$\dfrac{d\Phi}{dk}=0$，即

$$
r+2\frac{d\delta_0}{dk}-\frac{d\omega}{dk}\,t=0.
$$

引入群速度 $v_g=\dfrac{d\omega}{dk}$，得

$$
t=\frac{r+2\dfrac{d\delta_0}{dk}}{v_g}.
$$

### 3.3 时间延迟

若无相互作用（$\delta_0=0$），同样位置的波包中心到达时刻为 $t_0=\dfrac{r}{v_g}$。二者之差即**时间延迟**：

$$
\Delta t=t-t_0=\frac{2}{v_g}\frac{d\delta_0}{dk}.
$$

用 $\dfrac{1}{v_g}\dfrac{d}{dk}=\dfrac{d}{d\omega}=\hbar\dfrac{d}{dE}$，得

$$
\boxed{\;\Delta t=2\hbar\,\frac{d\delta_0}{dE}\;}
$$

这就是 **Eisenbud–Wigner 时间延迟公式**。因子 $2$ 来自 $S_0=e^{2i\delta_0}$ 中的 $2\delta_0$。

### 3.4 物理意义

- $\Delta t>0$：粒子"滞留"在散射区（如共振、准束缚态），出射延迟。
- $\Delta t<0$：粒子"超前"出射（如纯排斥势把粒子"推出"）。

但时间超前不能任意大——这正是下一节因果性的内容。

---

## 4. 因果性原理

### 4.1 物理陈述

因果性要求：**散射波不能在入射波到达散射体之前就离开散射体**。换言之，时间超前有一个由相互作用程 $R$ 决定的下界。粗略地，最大时间超前约为穿越相互作用区的时间：

$$
\Delta t\gtrsim-\frac{2R}{v_g},
$$

即入射波从 $r=-R$ 进入、再从 $r=+R$ 出去，至少需要时间 $2R/v_g$。散射波不可能比这更早出现。

### 4.2 与解析性的关系

因果性等价于散射振幅作为能量 $E$ 的函数在**上半平面解析**（无超前意味着推迟的格林函数只含正频率传播）。这一解析性又通过 Kramers–Kronig 关系把相移的实部与虚部（吸收）联系起来。本文不展开，但这是因果性的另一面。

### 4.3 因果性的量子力学精确表述

Wigner 发现，因果性可以严格表述为：径向波函数的**倒数对数导数**在固定半径 $r=R$ 处随能量**单调不减**。下一节从薛定谔方程严格证明这一点。

---

## 5. 倒数对数导数的单调性（Wronskian 严格推导）

### 5.1 两个不同能量的解

设 $u_1,u_2$ 分别是能量 $E_1,E_2$ 的两个径向解，都满足 $u(0)=0$：

$$
u_1''-\frac{2m}{\hbar^2}\big[V(r)-E_1\big]u_1=0,
\qquad
u_2''-\frac{2m}{\hbar^2}\big[V(r)-E_2\big]u_2=0.
$$

### 5.2 Wronskian 恒等式

用 $u_2$ 乘第一个方程、$u_1$ 乘第二个方程，相减：

$$
u_2u_1''-u_1u_2''=\frac{2m}{\hbar^2}(E_2-E_1)\,u_1u_2 .
$$

左边恰是导数 $\dfrac{d}{dr}\big[u_2u_1'-u_1u_2'\big]$：

$$
\frac{d}{dr}\big[u_2u_1'-u_1u_2'\big]=\frac{2m}{\hbar^2}(E_2-E_1)\,u_1u_2 .
$$

从 $0$ 到 $R$ 积分。在 $r=0$ 处 $u_1(0)=u_2(0)=0$，故

$$
u_2(R)u_1'(R)-u_1(R)u_2'(R)=\frac{2m}{\hbar^2}(E_2-E_1)\int_0^R u_1u_2\,dr .
$$

### 5.3 对数导数的能量导数

两端除以 $u_1(R)u_2(R)$：

$$
\frac{u_1'(R)}{u_1(R)}-\frac{u_2'(R)}{u_2(R)}
=\frac{2m}{\hbar^2}(E_2-E_1)\frac{\displaystyle\int_0^R u_1u_2\,dr}{u_1(R)u_2(R)} .
$$

定义**对数导数** $L(E)\equiv\dfrac{u'(R)}{u(R)}$（它只依赖能量，不依赖归一化）。则

$$
L(E_1)-L(E_2)=\frac{2m}{\hbar^2}(E_2-E_1)\frac{\displaystyle\int_0^R u_1u_2\,dr}{u_1(R)u_2(R)} .
$$

令 $E_1\to E_2=E$，$L(E_1)-L(E)\to(E_1-E)\dfrac{dL}{dE}$，右端 $\to\dfrac{2m}{\hbar^2}(E-E_1)\dfrac{\int_0^R u^2 dr}{u^2(R)}$。比较得

$$
\boxed{\;\frac{dL}{dE}=-\frac{2m}{\hbar^2}\frac{\displaystyle\int_0^R [u(r)]^2\,dr}{[u(R)]^2}\le 0\;}.
$$

即**对数导数随能量单调不增**。等价地，**倒数对数导数** $\mathcal R(E)\equiv\dfrac{u(R)}{u'(R)}=\dfrac1{L(E)}$ 满足

$$
\boxed{\;\frac{d\mathcal R}{dE}=-\frac{1}{L^2}\frac{dL}{dE}\ge 0\;}
$$

即**倒数对数导数随能量单调不减**。

### 5.4 物理解释

$\mathcal R(E)=u/u'$ 可理解为：把内部解 $u$ 线性外推到零（$u+u'(r-R)=0$）得到的"等效外推点"到原点的距离，即波函数节点的位置。随能量升高，波函数振荡更快，节点向内移动，$\mathcal R$ 减小——其单调性正是"波函数节点随能量单调移动"的数学表述。

这一单调性就是因果性：若 $\mathcal R$ 非单调，则在某个能量区间内，高能分量会"跑得比低能分量还慢"甚至反向，导致波包出现违背因果的超前传播。

---

## 6. Wigner 因果界的推导

### 6.1 匹配条件

对 $r>R$，$V=0$，解为自由波。$s$ 波解（归一化 $u=\dfrac{\sin(kr+\delta_0)}{\sin\delta_0}$，使零能极限下回到 $1-\dfrac ra$）满足

$$
u(r)=\frac{\sin(kr+\delta_0)}{\sin\delta_0},\qquad
u'(r)=k\,\frac{\cos(kr+\delta_0)}{\sin\delta_0}.
$$

在 $r=R$ 处，内部解的倒数对数导数必须等于外部自由波的倒数对数导数（$u,u'$ 连续）：

$$
\frac{u(R)}{u'(R)}=\frac{\sin(kR+\delta_0)}{k\cos(kR+\delta_0)}
=\frac1k\tan(kR+\delta_0).
$$

### 6.2 施加因果性（单调性）

由第 5 节，$\dfrac{u(R)}{u'(R)}$ 随 $E$ 单调不减：

$$
\frac{d}{dE}\left[\frac1k\tan(kR+\delta_0)\right]\ge 0 .
$$

利用 $\dfrac{d}{dE}=\dfrac{m}{\hbar^2k}\dfrac{d}{dk}$，并记 $\phi\equiv kR+\delta_0$：

$$
\frac{m}{\hbar^2k}\left[-\frac1{k^2}\tan\phi+\frac1k\left(R+\frac{d\delta_0}{dk}\right)\sec^2\phi\right]\ge 0 .
$$

方括号内非负：

$$
-\frac1{k^2}\tan\phi+\frac1k\left(R+\frac{d\delta_0}{dk}\right)\sec^2\phi\ge 0 .
$$

乘以 $k^2$（$k>0$）：

$$
-\tan\phi+k\left(R+\frac{d\delta_0}{dk}\right)\sec^2\phi\ge 0 .
$$

两边除以 $\sec^2\phi$（$>0$）：

$$
k\left(R+\frac{d\delta_0}{dk}\right)\ge\frac{\tan\phi}{\sec^2\phi}=\sin\phi\cos\phi=\frac12\sin(2\phi).
$$

于是得到 **Wigner 因果界（相移形式）**：

$$
\boxed{\;\frac{d\delta_0}{dk}\ge -R+\frac{1}{2k}\sin\!\big[2(kR+\delta_0)\big]\;}
$$

这就是因果性对相移导数施加的精确下界。注意：等号右侧 $-R$ 对应最大时间超前 $\Delta t\approx-\dfrac{2R}{v_g}$（见第 4 节），与直观一致；$\sin$ 项是 $s$ 波干涉的精确修正。

---

## 7. 低能极限与有效程的 Wigner 界

### 7.1 有效程展开

低能 $s$ 波相移满足有效程展开：

$$
k\cot\delta_0=-\frac1a+\frac12 r_0k^2+O(k^4),
$$

其中 $a$ 为散射长度，$r_0$ 为有效程。反解得 $\cot\delta_0=-\dfrac{1}{ak}+\dfrac12r_0k+O(k^3)$，故

$$
\tan\delta_0=\frac{1}{\cot\delta_0}
=\frac{-ak}{1-\dfrac{ar_0k^2}{2}}
=-ak-\frac{a^2r_0}{2}k^3+O(k^5).
$$

### 7.2 相移的小 $k$ 展开（含 arctan 修正）

注意 $\delta_0=\arctan(\tan\delta_0)$，不能直接令 $\delta_0=\tan\delta_0$（差一个 $k^3$ 阶项）。展开 $\arctan x=x-\dfrac{x^3}{3}+O(x^5)$，$x=\tan\delta_0$：

$$
\tan^3\delta_0=(-ak)^3+O(k^5)=-a^3k^3+O(k^5),
$$

故

$$
\delta_0=\tan\delta_0-\frac{\tan^3\delta_0}{3}+O(k^5)
=-ak-\frac{a^2r_0}{2}k^3+\frac{a^3}{3}k^3+O(k^5),
$$

即

$$
\delta_0=-ak+a^2\left(\frac a3-\frac{r_0}{2}\right)k^3+O(k^5),
$$

$$
\frac{d\delta_0}{dk}=-a+a^2\left(a-\frac{3r_0}{2}\right)k^2+O(k^4).
$$

### 7.3 代入因果界并取零能极限

把展开代入第 6.2 节的界。左端：

$$
\frac{d\delta_0}{dk}=-a+a^2\left(a-\frac{3r_0}{2}\right)k^2+O(k^4).
$$

右端：记 $c=\dfrac a3-\dfrac{r_0}2$，则 $2(kR+\delta_0)=2k(R-a)+2a^2ck^3$。展开 $\sin x=x-\dfrac{x^3}{6}+\cdots$：

$$
\sin\!\big[2(kR+\delta_0)\big]
=2k(R-a)+2a^2ck^3-\frac{4(R-a)^3}{3}k^3+O(k^5),
$$

$$
\frac{1}{2k}\sin\!\big[2(kR+\delta_0)\big]
=(R-a)+a^2ck^2-\frac{2(R-a)^3}{3}k^2+O(k^4),
$$

于是右端

$$
-R+\frac{1}{2k}\sin\!\big[2(kR+\delta_0)\big]
=-a+a^2\left(\frac a3-\frac{r_0}{2}\right)k^2-\frac{2(R-a)^3}{3}k^2+O(k^4).
$$

因果界要求（比较 $k^2$ 阶系数，$k^0$ 阶两边相等为 $-a$）：

$$
a^2\left(a-\frac{3r_0}{2}\right)\ge a^2\left(\frac a3-\frac{r_0}{2}\right)-\frac{2(R-a)^3}{3}.
$$

两端除以 $a^2$ 并化简：

$$
a-\frac{3r_0}{2}\ge\frac a3-\frac{r_0}{2}-\frac{2(R-a)^3}{3a^2}
$$

$$
\frac{2a}{3}-r_0\ge-\frac{2(R-a)^3}{3a^2}
$$

$$
r_0\le\frac{2a}{3}+\frac{2(R-a)^3}{3a^2}
=\frac{2[a^3+(R-a)^3]}{3a^2}.
$$

展开 $(R-a)^3=R^3-3R^2a+3Ra^2-a^3$，$a^3+(R-a)^3=R^3-3R^2a+3Ra^2$，最终得到

$$
\boxed{\;r_0\le 2R-\frac{2R^2}{a}+\frac{2R^3}{3a^2}\;}
$$

这就是**有效程的 Wigner 因果界**（三维 $s$ 波）。

> **硬球验证**：对半径 $R$ 的硬球，$a=R$，$r_0=\dfrac{2R}{3}$（可由 Bethe 公式或相移直接算出）。代入界：右边 $=2R-\dfrac{2R^2}{R}+\dfrac{2R^3}{3R^2}=2R-2R+\dfrac{2R}{3}=\dfrac{2R}{3}$。界被饱和！硬球正是"最紧凑"的相互作用，恰好达到因果性允许的极限。∎

---

## 8. 用 Bethe 公式交叉验证

Wigner 界也可以从 Bethe 有效程公式独立导出，两条路线殊途同归。

### 8.1 Bethe 公式

从零能与有限能两个径向方程出发，用 Wronskian 技巧可得（推导见《Causality and universality》一文）

$$
r_0=2\int_0^\infty\left[\left(1-\frac ra\right)^2-u_0^2(r)\right]dr,
$$

其中 $u_0(r)$ 为零能 $s$ 波径向波函数，渐近归一 $u_0(r)\to1-\dfrac ra$。

### 8.2 由 Bethe 公式得界

由于 $r>R$ 时 $u_0(r)=1-\dfrac ra$（精确成立），被积函数在 $r>R$ 恒为零，故

$$
r_0=2\int_0^R\left[\left(1-\frac ra\right)^2-u_0^2(r)\right]dr
=2\int_0^R\left(1-\frac ra\right)^2dr-2\int_0^R u_0^2(r)\,dr.
$$

第一项可直接积出

$$
2\int_0^R\left(1-\frac ra\right)^2dr=2R-\frac{2R^2}{a}+\frac{2R^3}{3a^2}.
$$

第二项 $\displaystyle 2\int_0^R u_0^2 dr\ge 0$（半正定），故

$$
r_0\le 2R-\frac{2R^2}{a}+\frac{2R^3}{3a^2}.
$$

与第 7 节结果完全一致。这条路线揭示了 Wigner 界与"因果性"的深层联系：Bethe 公式中的非负积分项 $2\int_0^R u_0^2$ 正是被"扣除"掉的因果部分，扣除后有效程不可能超过 $b(R)=2R-\dfrac{2R^2}{a}+\dfrac{2R^3}{3a^2}$。

---

## 9. 推广到任意维度 $d$ 与角动量 $L$

Wigner 界可推广到任意维度与分波（详见 Hammer & Lee 2009）。在 $d$ 维、角动量 $L$ 时，径向方程（$\hbar=1$）为

$$
p^2u_{L,d}^{(p)}=\left[-\frac{d^2}{dr^2}+\frac{(2L+d-3)(2L+d-1)}{4r^2}\right]u_{L,d}^{(p)}
+2\mu\int_0^R dr'\,W(r,r')u_{L,d}^{(p)}(r'),
$$

其中 $p^2/(2\mu)=E$。有效程 $r_{L,d}$ 满足推广的 Bethe 公式

$$
r_{L,d}=b_{L,d}(r)-2\int_0^r\big[u^{(0)}_{L,d}(r')\big]^2dr'\qquad(r\ge R),
$$

由被积函数半正定得推广的 Wigner 界

$$
r_{L,d}\le b_{L,d}(r)\qquad(\text{任意 }r\ge R),
$$

其中 $b_{L,d}(r)$ 的一般形式为

$$
b_{L,d}(r)=-\frac{\Gamma(L+\tfrac d2-1)\Gamma(L+\tfrac d2-\tfrac12)}{\pi}\left(\frac r2\right)^{-2L-d+4}
-\frac{4}{L+\frac d2-1}\frac{1}{a_{L,d}}\left(\frac r2\right)^2
+\frac{2\pi}{\Gamma(L+\frac d2)\Gamma(L+\frac d2+1)}\frac{1}{a_{L,d}^2}\left(\frac r2\right)^{2L+d}.
$$

$L=0,d=3$ 时上式第一项（$\propto r^{-1}$）的正规处理回到 $b(r)=2r-\dfrac{2r^2}{a}+\dfrac{2r^3}{3a^2}$，即本文的三维 $s$ 波结果。

---

## 10. 物理意义与应用

1. **有效程不能任意为负**：因果性禁止把 $r_0$ 调到界以下。这解释了为何物理散射的有效程有一个由散射长度 $a$ 与程 $R$ 决定的下界；想取零程极限 $R\to0$ 就必须让 $r_0$ 足够负，而这是被禁止的。

2. **普适性的限度**：在三维 $s$ 波（$2L+d=3$），零程极限仍可只用散射长度 $a$ 描述（$r_0\to0$ 不被严格禁止）。但对 $2L+d\ge4$（如三维 $p$ 波 $2L+d=5$），Wigner 界在 $R\to0$ 时发散，**禁止把有效程调零**，于是必须引入第二个低能参数（散射体积 + 有效程）。这就是 Hammer & Lee 的核心结论，Wigner 界是它的数学基础。

3. **共振与时间延迟**：Wigner 界把可观测的散射相移（时间延迟）与相互作用范围联系起来，为从相移数据反推相互作用的最小程提供了工具（如 $\alpha$–中子散射中各分波道的最小 $R$ 估计）。

4. **有效场论的指导**：在有效场论 / 有效程展开中，Wigner 界约束了低能常数的取值范围，指示哪些算符在领头阶就必须保留，保证低能展开的自洽与因果。

5. **深层含义**：因果性（散射波不先于入射波出现）在量子力学中不是外加假设，而是薛定谔方程（推迟边界条件）的必然结果，最终表现为波函数对数导数单调、相移导数有下界这一可计算的约束。它是少有的把"时间箭头 / 因果"转化为定量不等式的例子。

---

**一句话总结**：散射时间延迟 $\Delta t=2\hbar\dfrac{d\delta}{dE}$ 不能任意为负（散射波不能先于入射波出现），等价于倒数对数导数随能量单调不减；从薛定谔方程的 Wronskian 恒等式严格推出这一单调性，再配合有效程展开即得 Wigner 因果界 $r_0\le 2R-\dfrac{2R^2}{a}+\dfrac{2R^3}{3a^2}$。
