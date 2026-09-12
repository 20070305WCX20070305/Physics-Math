# 低能量子散射中的因果性与普适性

> **原文**：H.-W. Hammer and Dean Lee, *Causality and universality in low-energy quantum scattering*, Phys. Lett. B **681**, 500 (2009); arXiv:0907.1763。

---

## 目录

1. 文章大意
2. 基本设定：$d$ 维径向薛定谔方程
3. 渐近形式与相移
4. 有效程展开
5. Wronskian 恒等式
6. Bethe 积分公式的推广
7. Wigner 因果界
8. 低能普适性：一参数还是两参数？
9. 浅束缚态的波函数尺寸
10. 范德瓦尔斯尾的修正
11. 意义

---

## 1. 文章大意

本文把 **Wigner 因果界** 与 **Bethe 有效程积分公式** 推广到任意维度 $d$、任意角动量 $L$。

核心结论：

1. 因果性要求散射波不会在入射波到达散射体之前就传播出来，这给有效程 $r_{L,d}$ 一个上界（Wigner 界）：

   $$
   r_{L,d}\le b_{L,d}(r)\qquad(\text{对任意 } r\ge R).
   $$
2. 当 $2L+d\le3$ 时，低能散射只需要**一个**低能参数（三维 $s$ 波只需散射长度 $a$）。
   当 $2L+d\ge4$ 时，因果性禁止把有效程调到零，必须引入**第二个**低能参数（如三维 $p$ 波需要散射体积与有效程两个参数）。
3. 对浅束缚态，这第二个参数还刻画束缚态波函数的空间尺寸。
4. 对碱金属原子的范德瓦尔斯尾，需对有效程公式做相应修正。

---

## 2. 基本设定：$d$ 维径向薛定谔方程

### 2.1 从薛定谔方程到径向方程

考虑 $d$ 维空间中两个无自旋粒子，约化质量 $\mu$，转动不变的两体相互作用。在质心系中取 $\hbar=1$，能量 $E=p^2/(2\mu)$。

$d$ 维拉普拉斯算子作用在径向函数 $R_{L,d}(r)$ 上：

$$
\nabla^2 R=\frac{1}{r^{d-1}}\frac{d}{dr}\left(r^{d-1}\frac{dR}{dr}\right)-\frac{L(L+d-2)}{r^2}R,
$$

其中 $L$ 为 $d$ 维角动量量子数（离心势 $L(L+d-2)/r^2$）。

### 2.2 重新标度径向波函数

定义

$$
\boxed{\;u_{L,d}^{(p)}(r)=(pr)^{(d-1)/2}R_{L,d}^{(p)}(r)\;}
$$

注意这里把 $(pr)^{(d-1)/2}$ 写作 $p^{(d-1)/2}r^{(d-1)/2}$，在固定的 $p$ 下相当于 $r^{(d-1)/2}$（只差一个与 $r$ 无关的常数 $p^{(d-1)/2}$，它的作用是使 $u$ 在 $r\ge R$ 处的渐近式与贝塞尔函数对应整齐，见第 3 节）。

### 2.3 消去一阶导数项

令 $u=r^{(d-1)/2}R$（暂时略去常数 $p^{(d-1)/2}$）。逐项计算：

$$
\frac{dR}{dr}=\frac{d}{dr}\big[r^{-(d-1)/2}u\big]=r^{-(d-1)/2}u'-\frac{d-1}{2}r^{-(d+1)/2}u,
$$

$$
r^{d-1}\frac{dR}{dr}=r^{(d-1)/2}u'-\frac{d-1}{2}r^{(d-3)/2}u,
$$

再求导：

$$
\frac{d}{dr}\left[r^{d-1}\frac{dR}{dr}\right]
=\underbrace{r^{(d-1)/2}u''+\frac{d-1}{2}r^{(d-3)/2}u'}_{\frac d{dr}[r^{(d-1)/2}u']}
-\underbrace{\left[\frac{d-1}{2}r^{(d-3)/2}u'+\frac{(d-1)(d-3)}{4}r^{(d-5)/2}u\right]}_{\frac{d-1}{2}\frac d{dr}[r^{(d-3)/2}u]}.
$$

两个 $u'$ 项相消，得

$$
\frac{1}{r^{d-1}}\frac{d}{dr}\left[r^{d-1}\frac{dR}{dr}\right]
=r^{-(d-1)/2}u''-\frac{(d-1)(d-3)}{4}r^{-(d+3)/2}u.
$$

代入径向方程，并同乘 $r^{(d-1)/2}$：

$$
-u''+\frac{(d-1)(d-3)}{4r^2}u+\frac{L(L+d-2)}{r^2}u=p^2 u .
$$

### 2.4 合并离心项

$$
\frac{(d-1)(d-3)}{4}+L(L+d-2)
=L^2+Ld-2L+\frac{d^2-4d+3}{4}.
$$

另一方面展开

$$
\frac{(2L+d-3)(2L+d-1)}{4}
=\frac{(2L+d-2)^2-1}{4}
=L^2+Ld-2L+\frac{d^2-4d+3}{4}.
$$

二者相等，故

$$
\boxed{\;-u''+\frac{(2L+d-3)(2L+d-1)}{4r^2}u=p^2u\;}
$$

即自由径向方程。加入相互作用（写成实对称算符，核 $W(r,r')$，有限程 $R$）：

$$
\boxed{\;p^2u_{L,d}^{(p)}(r)=\left[-\frac{d^2}{dr^2}+\frac{(2L+d-3)(2L+d-1)}{4r^2}\right]u_{L,d}^{(p)}(r)
+2\mu\int_0^R dr'\,W(r,r')u_{L,d}^{(p)}(r')\;}
$$

这是论文式 (2)。

> 短程正则性条件：要求 $u^{(p)}_{L,d}$ 在 $r\to0$ 时为零、$\dfrac{d}{dr}u^{(p)}_{L,d}$ 有限。对局域势 $W(r,r')=V(r)\delta(r-r')$，只要 $V(r)=O(r^{-2+\epsilon})$（$\epsilon>0$）即可。

---

## 3. 渐近形式与相移

### 3.1 自由方程的贝塞尔解

令

$$
\nu\equiv L+\frac d2-1 .
$$

则离心系数

$$
\frac{(2L+d-3)(2L+d-1)}{4}=\frac{(2\nu-1)(2\nu+1)}{4}=\nu^2-\frac14.
$$

自由方程（$r\ge R$ 时 $W=0$）为

$$
u''+\left(p^2-\frac{\nu^2-1/4}{r^2}\right)u=0.
$$

这正是 $\nu$ 阶贝塞尔方程。其解为 $\sqrt{r}J_\nu(pr)$ 与 $\sqrt{r}Y_\nu(pr)$（$J_\nu$ 为第一类、$Y_\nu$ 为第二类贝塞尔函数）。一般解：

$$
u^{(p)}(r)=\sqrt r\big[A J_\nu(pr)+B\,Y_\nu(pr)\big]\qquad(r\ge R).
$$

### 3.2 相移的定义与归一化

论文选择归一化（式 3）：

$$
\boxed{\;u^{(p)}_{L,d}(r)=\sqrt{\frac{rp\pi}{2}}\,p^{L+d/2-1/2}
\left[\cot\delta_{L,d}(p)\,J_{L+d/2-1}(pr)-Y_{L+d/2-1}(pr)\right]\;}
$$

即取 $A\propto\cot\delta_{L,d}$、$B\propto -1$，比例系数为 $\sqrt{\frac{\pi}{2}}p^{L+d/2}$（注意 $\sqrt{\frac{rp\pi}{2}}\,p^{L+d/2-1/2}=\sqrt{\frac{r\pi}{2}}\,p^{L+d/2}$）。这种"驻波 + 相移"的形式是标准约定：相移 $\delta_{L,d}$ 通过 $\cot\delta_{L,d}$ 进入解。

### 3.3 散射振幅

分波散射振幅与相移的关系为（式 4）

$$
\boxed{\;f_{L,d}(p)\propto\frac{p^{2L}}{p^{2L+d-2}\cot\delta_{L,d}(p)-i\,p^{2L+d-2}}\;}
$$

分母出现组合 $p^{2L+d-2}\cot\delta_{L,d}$，这正是有效程展开的对象。

---

## 4. 有效程展开

对有限程相互作用，$p^{2L+d-2}\cot\delta_{L,d}(p)$ 是 $p$ 的偶函数（对于奇 $d$ 是解析偶函数；对偶 $d$ 会出现对数项）。一般形式为（式 5）

$$
\boxed{\;p^{2L+d-2}\left[\cot\delta_{L,d}(p)-\delta_{(d\bmod 2),0}\,\frac{2}{\pi}\ln(p\rho_{L,d})\right]
=-\frac{1}{a_{L,d}}+\frac12 r_{L,d}p^2+\sum_{n=0}^{\infty}(-1)^{n+1}P^{(n)}_{L,d}p^{2n+4}\;}
$$

其中：

- $\delta_{(d\bmod2),0}$ 在 $d$ 偶时为 $1$、$d$ 奇时为 $0$（偶维度的对数项）。
- $a_{L,d}$ 为**散射参数**（三维 $s$ 波即散射长度）。
- $r_{L,d}$ 为**有效程参数**。
- $P^{(n)}_{L,d}$ 为第 $n$ 阶形状参数。
- $\rho_{L,d}$ 为任意长度尺度，重新标度它只移动 $p^{2L+d-2}$ 的常数项；记 $\bar\rho_{L,d}$ 为使该常数为零的特殊取值。

低能（$p\to0$）主项：

$$
p^{2L+d-2}\cot\delta_{L,d}\to-\frac{1}{a_{L,d}}.
$$

---

## 5. Wronskian 恒等式

设 $u^{(p)}$、$u^{(p')}$ 是动量为 $p$、$p'$ 的两个解，都满足式 (2)。写出两个方程：

$$
u^{(p)''}=V_c\,u^{(p)}+2\mu\,\hat W u^{(p)}-p^2u^{(p)},\qquad
u^{(p')''}=V_c\,u^{(p')}+2\mu\,\hat W u^{(p')}-p'^2u^{(p')},
$$

其中 $V_c=\dfrac{(2L+d-3)(2L+d-1)}{4r^2}$，$\hat W u=\int_0^R W(r,r')u(r')dr'$。

构造 Wronskian（论文式 6）

$$
\boxed{\;u^{(p)}_{L,d}\frac d{dr}u^{(p')}_{L,d}-u^{(p')}_{L,d}\frac d{dr}u^{(p)}_{L,d}\;},
$$

对其求导：

$$
\frac{d}{dr}\big[u^{(p)}u^{(p')'}-u^{(p')}u^{(p)'}\big]
=u^{(p)}u^{(p')''}-u^{(p')}u^{(p)''}
$$

$$
=\big(p^2-p'^2\big)u^{(p)}u^{(p')}+2\mu\big[u^{(p)}\hat W u^{(p')}-u^{(p')}\hat W u^{(p)}\big].
$$

从 $0$ 到 $r$ 积分。在 $r=0$ 处 $u^{(p)}(0)=u^{(p')}(0)=0$，故下限为零。对于 $r\ge R$，相互作用项积分为零：因 $\hat W u$ 只在 $r<R$ 内非零，且 $W$ 为实对称核，故

$$
\int_0^r\big[u^{(p)}\hat W u^{(p')}-u^{(p')}\hat W u^{(p)}\big]dr'
=\langle u^{(p)}|W|u^{(p')}\rangle-\langle u^{(p')}|W|u^{(p)}\rangle=0 .
$$

于是得到关键的 Wronskian 恒等式：

$$
\boxed{\;u^{(p)}(r)u^{(p')'}(r)-u^{(p')}(r)u^{(p)'}(r)
=(p^2-p'^2)\int_0^r u^{(p)}(r')u^{(p')}(r')dr'\qquad(r\ge R)\;}
$$

---

## 6. Bethe 积分公式的推广

### 6.1 三维 $s$ 波：恢复 Bethe 原始公式

先看 $d=3,L=0$（$\nu=1/2$）。零能解 $u^{(0)}(r)$ 与有限能解 $u^{(p)}(r)$ 都满足自由方程（$r\ge R$），且由式 (3) 归一化：

- 零能（$p\to0$，$p\cot\delta\to-1/a$）解渐近：
  $$
  u^{(0)}(r)=1-\frac{r}{a},
  $$
- 有限能解渐近：
  $$
  u^{(p)}(r)=\frac{\sin(pr+\delta_0)}{\sin\delta_0}.
  $$

把 $p'=0$ 代入 Wronskian 恒等式：

$$
u^{(p)}(r)u^{(0)'}(r)-u^{(0)}(r)u^{(p)'}(r)=p^2\int_0^r u^{(p)}u^{(0)}dr'.
$$

左端用渐近式：$u^{(0)'}=-1/a$，$u^{(p)}=\dfrac{\sin(pr+\delta_0)}{\sin\delta_0}$，$u^{(p)'}=p\dfrac{\cos(pr+\delta_0)}{\sin\delta_0}$：

$$
\begin{aligned}
\text{左端}&=-\frac1a\frac{\sin(pr+\delta_0)}{\sin\delta_0}-p\frac{\cos(pr+\delta_0)}{\sin\delta_0}\left(1-\frac ra\right).
\end{aligned}
$$

用 $\cot\delta_0$ 展开：记 $\sin(pr+\delta_0)=\sin(pr)\cos\delta_0+\cos(pr)\sin\delta_0$，$\cos(pr+\delta_0)=\cos(pr)\cos\delta_0-\sin(pr)\sin\delta_0$，整理成

$$
\text{左端}=-\cot\delta_0\,A-B,
$$

其中

$$
A=\frac{\sin(pr)}{a}+p\left(1-\frac ra\right)\cos(pr),\qquad
B=\frac{\cos(pr)}{a}-p\left(1-\frac ra\right)\sin(pr).
$$

小 $p$ 展开（$\sin x=x-\frac{x^3}{6}+O(x^5)$，$\cos x=1-\frac{x^2}{2}+O(x^4)$）：

$$
A=p-p^3\left(\frac{r^2}{2}-\frac{r^3}{3a}\right)+O(p^5),\qquad
B=\frac1a-p^2\left(r-\frac{r^2}{2a}\right)+O(p^4).
$$

有效程展开 $p\cot\delta_0=-\dfrac1a+\dfrac12r_0p^2+O(p^4)$ 给出

$$
\cot\delta_0=-\frac1{pa}+\frac12r_0p+O(p^3).
$$

代入左端，$O(p^0)$ 项相消，$p^2$ 阶系数为

$$
\text{左端}=p^2\left[r-\frac{r^2}{a}+\frac{r^3}{3a^2}-\frac12r_0\right]+O(p^3).
$$

右端 $p^2\int_0^r u^{(p)}u^{(0)}dr'\to p^2\int_0^r\big[u^{(0)}\big]^2dr'$。比较 $p^2$ 系数：

$$
r-\frac{r^2}{a}+\frac{r^3}{3a^2}-\frac12r_0=\int_0^r\big[u^{(0)}(r')\big]^2dr'.
$$

左边恰是 $\displaystyle\int_0^r\left(1-\frac{r'}{a}\right)^2dr'$，于是

$$
\frac12 r_0=\int_0^r\left[\left(1-\frac{r'}{a}\right)^2-\big(u^{(0)}(r')\big)^2\right]dr'.
$$

取 $r\ge R$，$r'>R$ 时被积函数恒为零，上限可延拓到 $\infty$，得到 **Bethe 公式**：

$$
\boxed{\;r_0=2\int_0^\infty\left[\left(1-\frac ra\right)^2-\big(u^{(0)}(r)\big)^2\right]dr\;}.
$$

### 6.2 任意 $d,L$：式 (7) 与 $b_{L,d}(r)$

把上面的步骤逐字搬到任意 $d,L$。零能解 $u^{(0)}_{L,d}$（由式 (3) 的 $p\to0$ 极限、去掉整体 $p$ 因子后归一化）在 $r\ge R$ 处为

$$
u^{(0)}_{L,d}(r)=\sqrt{\frac{\pi r}{2}}\left[\frac{\Gamma(\nu)2^\nu}{\pi r^\nu}-\frac{1}{a_{L,d}}\frac{r^\nu}{2^\nu\Gamma(\nu+1)}\right],\qquad \nu=L+\frac d2-1.
$$

> 验证：$\nu=1/2$（三维 $s$ 波）时 $\Gamma(1/2)=\sqrt\pi$，$\Gamma(3/2)=\sqrt\pi/2$，上式化为
>
> $$
> u^{(0)}(r)=1-\frac ra,
> $$
>
> 与前面一致。∎

Bethe 公式推广为

$$
\boxed{\;r_{L,d}=b_{L,d}(r)-2\int_0^r\left[u^{(0)}_{L,d}(r')\right]^2dr'\qquad(r\ge R)\;}
$$

其中 $b_{L,d}(r)$ 是"自由零能解平方的积分"，由贝塞尔函数展开求得。分三种情形：

**(i) $2L+d=2$**（如二维 $s$ 波，$\nu=0$），含对数：

$$
b_{L,d}(r)=\frac{2r^2}{\pi}\left(\left[\ln\left(\frac r{2\rho_{L,d}}\right)+\gamma-\frac12+\frac{\pi}{2a_{L,d}}\right]^2+\frac14\right).
$$

**(ii) $2L+d=4$**（如三维 $p$ 波，$\nu=1$），含对数：

$$
b_{L,d}(r)=\frac{4}{\pi}\left[\ln\left(\frac r{2\rho_{L,d}}\right)+\gamma\right]-\frac{4}{a_{L,d}}\left(\frac r2\right)^2+\frac{\pi}{a_{L,d}^2}\left(\frac r2\right)^4.
$$

**(iii) 一般情形**（$2L+d$ 为任意正奇数或 $\ge6$ 的偶数）：

$$
\boxed{\;b_{L,d}(r)=-\frac{\Gamma(L+\tfrac d2-1)\Gamma(L+\tfrac d2-\tfrac12)}{\pi}\left(\frac r2\right)^{-2L-d+4}
-\frac{4}{L+\frac d2-1}\frac{1}{a_{L,d}}\left(\frac r2\right)^2
+\frac{2\pi}{\Gamma(L+\frac d2)\Gamma(L+\frac d2+1)}\frac{1}{a_{L,d}^2}\left(\frac r2\right)^{2L+d}\;}
$$

（$\gamma$ 为欧拉–马歇罗尼常数。）其中 $L=0,d=3$ 情形即 Bethe 原始公式（$b(r)=2r-\frac{2r^2}{a}+\frac{2r^3}{3a^2}$，见 6.1 的等价形式），$L$ 为一般值的推广由 Madsen 给出，这里推广到任意 $d$。

> **关键结构**：第一项 $\propto(r/2)^{-2L-d+4}$。当 $2L+d\ge5$ 时指数为负，且系数为负，故 $r\to0$ 时 $b_{L,d}(r)\to-\infty$。这是下一节 Wigner 界和普适性结论的根源。

---

## 7. Wigner 因果界

由式 (7)，被积函数 $[u^{(0)}]^2\ge0$ 半正定，故

$$
\boxed{\;r_{L,d}\le b_{L,d}(r)\qquad(\text{对任意 }r\ge R)\;}
$$

这是论文式 (11)。

**因果性解释**：散射时间延迟与相移的能量导数成正比

$$
\Delta t=2\hbar\frac{d\delta}{dE}.
$$

入射波必须先进入相互作用区域，散射波才能出来，故 $\dfrac{d\delta}{dE}$ 不能任意为负（否则出现任意大的时间超前）。因果性的量子力学表述是：倒数对数导数 $u^{(p)}/\dfrac d{dr}u^{(p)}$ 对能量导数非负。这一事实可由 Wronskian 恒等式导出。零能时 $\dfrac{d\delta_{L,d}}{dE}$ 正比于有效程 $r_{L,d}$，于是因果性给出 $r_{L,d}$ 的上界。

对 $\alpha$–中子散射的 $S_{1/2},P_{1/2},P_{3/2}$ 道，画出 $b_{L,3}(r)-r_{L,3}$，非负性要求给出 $R$ 的最小值：$S_{1/2}$ 约 $0.1$ fm，$P_{1/2}$ 约 $0.6$ fm，$P_{3/2}$ 约 $1.1$ fm（与 $\alpha$ 均方根半径、$\pi$ 介子康普顿波长约 $1.5$ fm 比较，说明有效场论截断尺度选取需谨慎）。

---

## 8. 低能普适性：一参数还是两参数？

固定相互作用程 $R$，取低能极限 $p\to0$：

- **$2L+d\le3$**：低能散射只依赖**一个**有量纲参数。

  - $2L+d=1$、$2L+d=3$：参数为 $a_{L,d}$（三维 $s$ 波即散射长度 $a$）。
  - $2L+d=2$：参数为 $\bar\rho_{L,d}$。
- **$2L+d\ge4$**：低能非微扰极限出现**第二个**有量纲参数。

在 $|a_{L,d}|\to\infty$ 极限下，式 (11) 的上界化为（式 12）：

$$
\bar\rho_{L,d}\le\frac r2 e^{\gamma}\quad(2L+d=4),\qquad
r_{L,d}\le-\frac{2\Gamma(L+\tfrac d2-1)\Gamma(L+\tfrac d2-\tfrac12)}{\pi}\left(\frac r2\right)^{-2L-d+4}\quad(2L+d\ge5).
$$

对 $2L+d=4$，无法通过微调参数压低 $p^2\ln(p\bar\rho_{L,d})$ 项，因为界禁止把对数宗量调到 1；对 $2L+d\ge5$，式 (12) 的负上界禁止把 $r_{L,d}$ 调为零以消去 $\frac12 r_{L,d}p^2$ 项。

因此在两种情形下，低能非微扰极限都剩下**两个相关参数**，对应重整化群不动点附近的两个相关方向。例如三维 $p$ 波中子–$\alpha$ 散射需在领头阶同时包含散射体积与有效程两个算符。

**一句话**：对于 $2L+d\ge4$，因果波传播对有限程相互作用达到标度不变极限构成**基本障碍**——有效程作为第二个相关低能参数必然涌现，无法在不违反因果性的前提下调零。

---

## 9. 浅束缚态的波函数尺寸

考虑 $p=i p_I$（虚动量）处的束缚态，零结合能极限 $p_I\to0^+$。记 $P_>(r)$ 为两粒子间距大于 $r$ 的概率：

$$
P_>(r)=\int_r^\infty dr'\left[\hat u^{(ip_I)}_{L,d}(r')\right]^2,
$$

其中 $\hat u^{(ip_I)}$ 为归一化束缚态波函数。

- 对 $2L+d\le3$：$p_I\to0^+$ 时 $P_>(r)\to1$（对任意 $r$）——束缚态波函数弥散到很大的距离，短程物理在低能下无关。
- 对 $2L+d\ge4$ 则不同：
  - $2L+d=4$：$P_>(r)$ 对数依赖 $\bar\rho_{L,d}$，可调到 $[0,1]$ 内任意值。
  - $2L+d\ge5$（式 14）：

$$
\boxed{\;P_>(r)\to\frac{\Gamma(L+\tfrac d2-1)\Gamma(L+\tfrac d2-\tfrac12)}{(-r_{L,d})\pi}\left(\frac r2\right)^{-2L-d+4}\qquad(r\ge R)\;}
$$

此时束缚态波函数的特征尺寸为

$$
\big(-r_{L,d}\big)^{1/(-2L-d+4)}.
$$

> **推导要点**：零结合能极限下，束缚态在 $r\ge R$ 处的形状趋于零能散射解 $u^{(0)}$（相差归一化），而这一归一化由 Wronskian 恒等式与有效程 $r_{L,d}$ 联系起来。故第二个低能参数 $r_{L,d}$ 同时刻画了浅束缚态的尺寸。

---

## 10. 范德瓦尔斯尾的修正

碱金属原子的 $p$ 波 Feshbach 共振需计入长程范德瓦尔斯相互作用

$$
W(r,r')=-C_6 r^{-6}\delta(r-r')\qquad(r,r'\ge R).
$$

引入长度尺度 $\beta_6=(2\mu C_6)^{1/4}$，并取 $d=3$。此时散射态应与 $r^{-6}$ 势的精确解比较，而非自由贝塞尔函数；$r<R$ 内的相互作用由解析于 $p^2$ 的有限程 $K$ 矩阵描述：

$$
K_L(p^2)=\sum_{n=0,1,\ldots}K_L^{(2n)}p^{2n}.
$$

相对自由球贝塞尔函数定义相移时，有效程展开不再解析于 $p^2$：$L=0$ 时领头非解析项 $\propto p^3$，$L=1$ 时 $\propto p^1$（使通常的有效程定义失效）。但对纯范德瓦尔斯尾，在零能共振极限（$K_L^{(0)}\to0$）下，$L=1$ 的 $p^1$ 系数也消失，于是 $s$ 波与 $p$ 波都可定义有效程（式 15）：

$$
r_0=\frac{[\Gamma(1/4)]^2}{3\pi}\big(\beta_6+3K_0^{(2)}\beta_6^{-3}\big),\qquad
r_1=-\frac{6[\Gamma(3/4)]^2}{5\pi}\big(\beta_6^{-3}-K_1^{(2)}\beta_6^{-1}\big).
$$

对碱金属原子单道散射，$K_L^{(2)}$ 相比 $\beta_6^2$ 可忽略；此时把 $R\sim\beta_6$ 代入式 (12) 的界，$L=0,1$ 均满足。若多道系统中 $K_L^{(2)}$ 不可忽略，则 $K_L^{(2)}$ 应满足类似的有效程 Wigner 界——这为研究多道 Feshbach 共振提供了出发点。

---

## 11. 意义

1. **普适性的完整判据**：本文第一次系统地给出"低能散射需要几个参数"的一般判据（$2L+d$ 与 $3$ 比较），统一了三维 $s$ 波（一参数）、三维 $p$ 波（两参数）等已知结果，并推广到任意维度与分波。
2. **因果性与重整化群的结合**：Wigner 因果界被翻译成重整化群语言——第二个相关方向（有效程）的涌现源于因果性对零程极限的禁止，把"为什么不能取零程"这一老问题讲清楚了。
3. **有效场论的指导**：对三维 $p$ 波散射（如中子–$\alpha$、$p$ 波 Feshbach 共振），指出领头阶就必须同时包含散射体积与有效程算符，纠正了只用散射体积的不足，直接指导有效场论/有效程展开的构造。
4. **束缚态几何的信息**：把有效程与浅束缚态的空间尺寸联系起来（式 14），为从散射可观测量提取束缚态结构提供了桥梁。
5. **长程力的推广**：范德瓦尔斯尾的分析表明有限程假设如何系统修正，连接了冷原子物理与核物理中的同类问题。

---

**一句话总结**：因果性要求散射波不在入射波之前出现，等价于有效程 $r_{L,d}$ 的上界（Wigner 界 $r_{L,d}\le b_{L,d}(r)$）；当 $2L+d\ge4$ 时该界在零程极限发散，从而禁止把有效程调零——低能散射必须由**两个**参数描述，且第二个参数同时决定浅束缚态的尺寸。
