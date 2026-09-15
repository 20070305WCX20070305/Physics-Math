# 横向约束下的原子散射与不可穿透玻色子气体（Tonks 气体）

> **原文**：M. Olshanii, *Atomic Scattering in the Presence of an External Confinement and a Gas of Impenetrable Bosons*, Phys. Rev. Lett. **81**, 938 (1998).

---

## 目录

1. 文章大意
2. 物理模型与基本薛定谔方程
3. 黄氏赝势与正规化部分
4. 横向模式展开与投影
5. Green 函数求解模式方程
6. 自洽条件与 $\mathcal{L}$ 函数的渐近展开
7. 一维散射振幅与一维散射长度
8. 低能极限：有效一维 $\delta$ 势
9. Tonks 气体与玻色–费米对偶
10. 意义

---

## 1. 文章大意

冷原子被限制在一个"雪茄形"波导或原子导波管中：横向（$X$-$Y$ 平面）由强简谐势阱约束，纵向（$Z$ 方向）自由。当横向运动冻结在基态、纵向动能低于横向激发能隙时，三维两体碰撞可以约化为一个**有效一维问题**。

文章的核心结论是：

1. 给出精确的一维偶宇称散射振幅
   $$
   f_{\rm even}(k_z)=-\frac{1}{1+ik_z a_{1D}-\dfrac{ik_z a_\perp}{2}\,\overline{\mathcal L}\!\left(-\dfrac{k_z^2a_\perp^2}{4}\right)}
   $$
   其中一维散射长度为
   $$
   a_{1D}=-\frac{a_\perp^2}{2a}\left(1-C\,\frac{a}{a_\perp}\right),\qquad C=1.4603\ldots=-\zeta(1/2).
   $$

2. 低能极限 $k_z|a_{1D}|\ll 1$ 下，散射退化为**全反射**，对应一维不可穿透玻色子（Tonks 气体）。

3. 相互作用可近似为有效一维 $\delta$ 势
   $$
   U_{1D}(z)=g_{1D}\,\delta(z),\qquad g_{1D}=-\frac{\hbar^2}{\mu a_{1D}}.
   $$

4. Tonks 气体与理想费米气体存在一一映射（玻色–费米对偶），动量分布出现 $1/\sqrt{k_z}$ 峰。

---

## 2. 物理模型与基本薛定谔方程

### 2.1 模型假设

记原子质量为 $m$，约化质量为
$$
\mu=\frac{m}{2}.
$$

- **(a)** 横向约束用轴对称二维简谐势近似，频率 $\omega_\perp$。
- **(b)** 沿 $Z$ 轴自由。
- **(c)** 原子间相互作用用黄氏赝势（见第 3 节）描述。
- **(d)** 原子运动（横向与纵向）被"冷却"到横向振动能 $\hbar\omega_\perp$ 以下。

### 2.2 质心与相对运动的分离

简谐约束允许分离质心与相对运动。相对坐标 $\mathbf r=\mathbf r_2-\mathbf r_1$ 满足的薛定谔方程为

$$
\boxed{\;\left[\frac{\hat p_z^2}{2\mu}+g\,\delta(\mathbf r)\,\frac{\partial}{\partial r}(r\,\cdot\,)+\hat H_\perp(\hat p_x,\hat p_y,x,y)\right]\Psi=E\Psi\;}
$$

其中横向哈密顿量为二维简谐振子

$$
\hat H_\perp=\frac{\hat p_x^2+\hat p_y^2}{2\mu}+\frac{\mu\omega_\perp^2}{2}(x^2+y^2),
$$

其本征谱为
$$
E_{n,m_z}=\hbar\omega_\perp(n+1),\qquad n=0,1,2,\ldots
$$

$n$ 为主量子数，$m_z$ 为绕 $Z$ 轴的角动量，$m_z=0,\pm2,\pm4,\ldots$（$n$ 偶）或 $m_z=\pm1,\pm3,\ldots$（$n$ 奇）。

二维谐振子本征态记为 $\phi_{n,m_z}(\rho,\varphi)$，其中 $\rho=\sqrt{x^2+y^2}$。基态波函数为

$$
\phi_{0,0}(\rho)=\frac{1}{\sqrt{\pi}\,a_\perp}\,e^{-\rho^2/2a_\perp^2},\qquad
a_\perp=\sqrt{\frac{\hbar}{\mu\omega_\perp}},
$$

于是
$$
\phi_{0,0}(0)=\frac{1}{\sqrt{\pi}\,a_\perp},\qquad |\phi_{0,0}(0)|^2=\frac{1}{\pi a_\perp^2}.
$$

一个关键性质（后面反复用到）：对**所有** $m_z=0$ 的态（即 $n$ 为偶数），都有

$$
\boxed{\;|\phi_{n,0}(0)|^2=\frac{1}{\pi a_\perp^2}\quad(\forall\ n=0,2,4,\ldots)\;}
$$

> **证明**：二维各向同性谐振子 $m_z=0$ 的径向波函数可写为
> $$
> \phi_{n,0}(\rho)=\frac{1}{\sqrt{2\pi}}\,R_{n,0}(\rho),\qquad
> R_{n,0}(\rho)=\sqrt{\frac{2}{a_\perp^2}}\,L_{n/2}\!\left(\frac{2\rho^2}{a_\perp^2}\right)e^{-\rho^2/2a_\perp^2},
> $$
> 其中 $L_k$ 为拉盖尔多项式，$L_k(0)=1$。故
> $$
> \phi_{n,0}(0)=\frac{1}{\sqrt{2\pi}}\sqrt{\frac{2}{a_\perp^2}}=\frac{1}{\sqrt{\pi}\,a_\perp},
> $$
> 与 $n$ 无关。∎

### 2.3 入射波与能量条件

设入射波处于横向基态
$$
\Psi_{\rm in}=e^{ik_z z}\,\phi_{0,0}(\rho).
$$

为保证碰撞后横向态不改变（即不发生向 $n>0$ 模式的跃迁），纵向动能需小于基态与第一个轴对称激发态（$n=2,m_z=0$）之间的能隙：

$$
\frac{\hbar^2 k_z^2}{2\mu}<E_{n=2,m_z=0}-E_{n=0,m_z=0}=2\hbar\omega_\perp .
$$

> 注意：碰撞过程中对高能轴对称模式（$n>0,m_z=0$）的**虚激发**并不被禁止，必须计入，这正是下面模式展开所做的事。

### 2.4 散射波函数的渐近形式

$$
\Psi(z,\rho)\xrightarrow{|z|\to\infty}
\Big\{e^{ik_z z}+f_{\rm even}\,e^{ik_z|z|}+f_{\rm odd}\,\mathrm{sgn}(z)\,e^{ik_z|z|}\Big\}\phi_{0,0}(\rho),
$$

$f_{\rm even}$ 与 $f_{\rm odd}$ 分别为一维偶、奇宇称散射振幅。本文将证明 $f_{\rm odd}=0$，并解析求出 $f_{\rm even}$。

---

## 3. 黄氏赝势与正规化部分

### 3.1 赝势的定义

三维相互作用用黄氏赝势近似：
$$
U(\mathbf r)=g\,\delta(\mathbf r)\,\frac{\partial}{\partial r}(r\,\cdot\,),\qquad
g=\frac{2\pi\hbar^2 a}{\mu},
$$

其中 $a$ 为"真实"三维相互作用势的 $s$ 波散射长度。

算子 $\dfrac{\partial}{\partial r}(r\,\cdot\,)$ 称为**正规化算子**，它把散射波中的 $1/r$ 发散去掉。

### 3.2 正规化部分 $\mathcal R$

赝势的作用效果是：若波函数在原点附近具有形式
$$
\Psi(\mathbf r)\xrightarrow{r\to0} A\left(\frac1r-\frac1a\right)+O(r),
$$
则
$$
\frac{\partial}{\partial r}\big[r\Psi(\mathbf r)\big]_{r\to0}=\frac{\partial}{\partial r}\big[A(1-r/a)\big]_{r\to0}=-\frac{A}{a}.
$$

我们定义**正规化部分**
$$
\boxed{\;\mathcal R\equiv\left.\frac{\partial}{\partial r}\big[r\Psi(\mathbf r)\big]\right|_{r\to0}
=\left.\frac{\partial}{\partial z}\big[z\Psi(z,\rho=0)\big]\right|_{z\to0^+}\;}
$$

（第二等式来自：沿 $Z$ 轴 $r=|z|$，偶宇称部分在 $z\to0^+$ 时 $r=z$。）

于是 $A=-a\mathcal R$，原点附近
$$
\Psi(\mathbf r)\sim -\frac{a\mathcal R}{r}+\mathcal R+\cdots,
$$

即奇异的 $1/r$ 项系数为 $-a\mathcal R$，正规部分（常数项）恰为 $\mathcal R$。

赝势作用于 $\Psi$ 给出
$$
U\Psi=g\,\delta(\mathbf r)\,\mathcal R.
$$

> **验证 $g=2\pi\hbar^2a/\mu$ 正确**：在无约束的自由三维空间中（$E=\hbar^2k^2/2\mu$），自由 Green 函数为
> $$
> G_0(\mathbf r)=-\frac{\mu}{2\pi\hbar^2}\frac{e^{ikr}}{r}.
> $$
> Lippmann–Schwinger 方程给出 $\Psi=\Psi_0+g\mathcal R\,G_0$，故
> $$
> \Psi\approx \Psi_0(0)-\frac{\mu g\mathcal R}{2\pi\hbar^2}\Big(\frac1r+ik\Big),\qquad (r\to0),
> $$
> 从而 $r\Psi\approx\Psi_0(0)r-\frac{\mu g\mathcal R}{2\pi\hbar^2}(1+ikr)$，
> $$
> \mathcal R=\partial_r(r\Psi)|_{r\to0}=\Psi_0(0)-\frac{i\mu g k}{2\pi\hbar^2}\mathcal R.
> $$
> 解得 $\mathcal R=\dfrac{\Psi_0(0)}{1+i\mu gk/2\pi\hbar^2}$。散射波 $\Psi_{\rm sc}=-\dfrac{\mu g\mathcal R}{2\pi\hbar^2}\dfrac{e^{ikr}}{r}$，故 $s$ 波振幅
> $$
> f_0=-\frac{\mu g\mathcal R}{2\pi\hbar^2}=-\frac{a}{1+iak}\Psi_0(0),
> $$
> 其中用了 $\dfrac{\mu g}{2\pi\hbar^2}=a$。这正是散射长度 $a$ 的 $s$ 波振幅。∎

---

## 4. 横向模式展开与投影

### 4.1 模式展开

由于入射波处于 $m_z=0$ 态，赝势是 $s$ 波（不传递角动量），故只需保留 $m_z=0$ 的横向模式（$n$ 为偶数）：

$$
\Psi(z,\rho)=\sum_{n=0,2,4,\ldots}\psi_n(z)\,\phi_{n,0}(\rho).
$$

### 4.2 投影到单个模式

把展开代入薛定谔方程。赝势项
$$
g\,\delta(\mathbf r)\,\mathcal R=g\,\delta(z)\,\delta^{(2)}(\boldsymbol\rho)\,\mathcal R,
$$

用 $\phi_{n,0}^*(\rho)$ 对横向坐标积分：
$$
\int d^2\rho\,\phi_{n,0}^*(\rho)\,g\,\delta(z)\delta^{(2)}(\boldsymbol\rho)\mathcal R
=g\,\delta(z)\,\mathcal R\,\phi_{n,0}^*(0).
$$

横向动能与势能给出 $E_{n,0}\psi_n$。总能量
$$
E=\frac{\hbar^2 k_z^2}{2\mu}+\hbar\omega_\perp
$$
（横向基态能 $\hbar\omega_\perp$ 加纵向动能）。于是第 $n$ 个模式满足

$$
\left[-\frac{\hbar^2}{2\mu}\frac{d^2}{dz^2}+E_{n,0}-E\right]\psi_n(z)=-g\,\mathcal R\,\phi_{n,0}^*(0)\,\delta(z).
$$

其中 $E_{n,0}-E=n\hbar\omega_\perp-\dfrac{\hbar^2k_z^2}{2\mu}$。

### 4.3 化为标准形式

两端乘 $-\dfrac{2\mu}{\hbar^2}$，并定义
$$
\kappa_n^2\equiv\frac{2\mu}{\hbar^2}\big(E_{n,0}-E\big)=\frac{2n}{a_\perp^2}-k_z^2,
$$

得到
$$
\boxed{\;\left[\frac{d^2}{dz^2}-\kappa_n^2\right]\psi_n^{({\rm sc})}(z)
=\frac{2\mu g}{\hbar^2}\,\mathcal R\,\phi_{n,0}^*(0)\,\delta(z)\;}
$$

这里 $\psi_n^{({\rm sc})}$ 是第 $n$ 个模式的散射部分（$n=0$ 模式还需加上入射波 $e^{ik_zz}$）。

---

## 5. Green 函数求解模式方程

### 5.1 Green 函数

方程 $\left[\dfrac{d^2}{dz^2}-\kappa_n^2\right]G_n(z)=\delta(z)$ 在"出射/衰减"边界条件下的解为

$$
G_n(z)=-\frac{e^{-\kappa_n|z|}}{2\kappa_n}.
$$

> **验证**：
> $$
> \frac{d}{dz}e^{-\kappa_n|z|}=-\kappa_n\,\mathrm{sgn}(z)\,e^{-\kappa_n|z|},
> \qquad
> \frac{d^2}{dz^2}e^{-\kappa_n|z|}=\kappa_n^2 e^{-\kappa_n|z|}-2\kappa_n\,\delta(z),
> $$
> 故
> $$
> \left[\frac{d^2}{dz^2}-\kappa_n^2\right]\left(-\frac{e^{-\kappa_n|z|}}{2\kappa_n}\right)
> =\delta(z).\ \square
> $$

### 5.2 各模式的散射波

$$
\psi_n^{({\rm sc})}(z)=\frac{2\mu g}{\hbar^2}\,\mathcal R\,\phi_{n,0}^*(0)
\left(-\frac{e^{-\kappa_n|z|}}{2\kappa_n}\right)
=-\frac{\mu g}{\hbar^2}\,\mathcal R\,\phi_{n,0}^*(0)\,\frac{e^{-\kappa_n|z|}}{\kappa_n}.
$$

对 $n=0$，$\kappa_0^2=-k_z^2$，取 $\kappa_0=-ik_z$（出射波 $e^{ik_z|z|}$），于是

$$
\psi_0^{({\rm sc})}(z)=-\frac{\mu g}{\hbar^2}\,\mathcal R\,\phi_{0,0}^*(0)\,\frac{e^{ik_z|z|}}{-ik_z}
=\frac{\mu g}{\hbar^2}\,\mathcal R\,\phi_{0,0}^*(0)\,\frac{e^{ik_z|z|}}{ik_z}.
$$

### 5.3 散射振幅与 $\mathcal R$ 的关系

由渐近定义 $\psi_0^{({\rm sc})}(z)=f_{\rm even}\,e^{ik_z|z|}$，得

$$
\boxed{\;f_{\rm even}=\frac{\mu g}{\hbar^2}\frac{\mathcal R\,\phi_{0,0}^*(0)}{ik_z}
=\frac{\gamma\,\mathcal R}{ik_z\,\phi_{0,0}(0)}\;}
$$

其中引进了无量纲（有量纲）组合
$$
\boxed{\;\gamma\equiv\frac{\mu g}{\hbar^2\pi a_\perp^2}=\frac{2a}{a_\perp^2}\;},
$$
（第二等式用 $g=2\pi\hbar^2a/\mu$ 得到。）

同时 $f_{\rm odd}=0$，因为赝势只产生 $s$ 波（偶宇称）散射。

### 5.4 波函数在 $\rho=0$ 处的表达式

$$
\Psi(z,0)=\phi_{0,0}(0)e^{ik_z z}+\sum_{n}\psi_n^{({\rm sc})}(z)\,\phi_{n,0}(0)
=\phi_{0,0}(0)e^{ik_z z}-\frac{\mu g\mathcal R}{\hbar^2}\sum_{n}|\phi_{n,0}(0)|^2\,\frac{e^{-\kappa_n|z|}}{\kappa_n}.
$$

利用 $|\phi_{n,0}(0)|^2=\dfrac{1}{\pi a_\perp^2}$，得

$$
\Psi(z,0)=\frac{1}{\sqrt\pi a_\perp}e^{ik_z z}
-\frac{\mu g\mathcal R}{\hbar^2\pi a_\perp^2}
\left[\frac{e^{ik_z|z|}}{-ik_z}+\sum_{n=2,4,\ldots}\frac{e^{-\kappa_n|z|}}{\kappa_n}\right].
$$

---

## 6. 自洽条件与 $\mathcal L$ 函数的渐近展开

### 6.1 分离 $n=0$ 模式并定义 $\mathcal L$

令 $n=2s'$（$s'=0,1,2,\ldots$）。对 $n\ge2$（$s'\ge1$）：

$$
\kappa_n^2=\frac{4s'}{a_\perp^2}-k_z^2=\frac{4}{a_\perp^2}\left(s'-\frac{k_z^2a_\perp^2}{4}\right)
=\frac{4}{a_\perp^2}(s'+\epsilon),
$$

其中
$$
\boxed{\;\epsilon\equiv-\frac{k_z^2a_\perp^2}{4}\;},\qquad
\kappa_n=\frac{2}{a_\perp}\sqrt{s'+\epsilon}.
$$

于是

$$
\sum_{n=2,4,\ldots}\frac{e^{-\kappa_n|z|}}{\kappa_n}
=\sum_{s'=1}^{\infty}\frac{\exp\!\left(-\frac{2|z|}{a_\perp}\sqrt{s'+\epsilon}\right)}
{\frac{2}{a_\perp}\sqrt{s'+\epsilon}}
=\frac{a_\perp}{2}\,\mathcal L\!\left(\frac{2|z|}{a_\perp},\epsilon\right),
$$

其中定义了

$$
\boxed{\;\mathcal L[\xi,\epsilon]\equiv\sum_{s'=1}^{\infty}\frac{e^{-\sqrt{s'+\epsilon}\,\xi}}{\sqrt{s'+\epsilon}}\;}
$$

于是

$$
\Psi(z,0)=\frac{1}{\sqrt\pi a_\perp}e^{ik_z z}
-\frac{\mu g\mathcal R}{\hbar^2\pi a_\perp^2}
\left[\frac{e^{ik_z|z|}}{-ik_z}+\frac{a_\perp}{2}\mathcal L\!\left(\frac{2|z|}{a_\perp},\epsilon\right)\right].
$$

### 6.2 $\mathcal L$ 的小 $\xi$ 展开

当 $\xi\to0$（即 $z\to0$）时，$\mathcal L$ 中的级数趋于发散（对应 $1/r$ 奇异）。用积分近似级数：

$$
\sum_{s'=1}^{\infty}\frac{e^{-\sqrt{s'+\epsilon}\,\xi}}{\sqrt{s'+\epsilon}}
\xrightarrow{\xi\to0}
\int_0^\infty\frac{e^{-\sqrt{s+\epsilon}\,\xi}}{\sqrt{s+\epsilon}}\,ds+\zeta\Big(\tfrac12\Big)+\cdots
$$

换元 $u=\sqrt{s+\epsilon}$，$ds=2u\,du$，$\dfrac{ds}{\sqrt{s+\epsilon}}=2du$：

$$
\int_0^\infty\frac{e^{-\sqrt{s+\epsilon}\,\xi}}{\sqrt{s+\epsilon}}ds
=\int_{\sqrt\epsilon}^{\infty}e^{-u\xi}\,2du
=\frac{2}{\xi}e^{-\sqrt\epsilon\,\xi}
=\frac{2}{\xi}+O(1).
$$

所以主导发散为 $\dfrac{2}{\xi}$。正规化后的常数部分：

$$
\mathcal L[\xi,\epsilon]=\frac{2}{\xi}-C+\overline{\mathcal L}(\epsilon)+O(\xi),
$$

其中

$$
\boxed{\;C=-\zeta\!\left(\tfrac12\right)=1.4603\ldots\;}
$$

其来源是
$$
C=\lim_{s\to\infty}\left[\int_0^s\frac{ds'}{\sqrt{s'}}-\sum_{s'=1}^{s}\frac{1}{\sqrt{s'}}\right],
$$
即积分 $\displaystyle\int_0^s\frac{ds'}{\sqrt{s'}}=2\sqrt s$ 与部分和 $\displaystyle\sum_{s'=1}^s\frac1{\sqrt{s'}}\approx 2\sqrt s+\zeta(1/2)$ 之差在 $s\to\infty$ 时的极限 $=-\zeta(1/2)$。

而 $\overline{\mathcal L}(\epsilon)$ 为 $\epsilon$ 的解析函数，$\overline{\mathcal L}(0)=0$，其级数形式为

$$
\overline{\mathcal L}(\epsilon)=\sum_{n=1}^{\infty}\frac{(-1)^n\,\zeta[(1+2n)/2]\,(2n-1)!!}{2^n\,n!}\,\epsilon^n,
$$

故 $\overline{\mathcal L}(\epsilon)=O(\epsilon)=O(k_z^2a_\perp^2)$。

### 6.3 自洽方程

由 $\mathcal R=\left.\partial_z\big[z\Psi(z,0)\big]\right|_{z\to0^+}$，代入 $\Psi(z,0)$ 的表达式逐项求导：

- 入射项 $z e^{ik_zz}$：$\partial_z(z e^{ik_zz})|_{z\to0}=1$，乘以 $\dfrac{1}{\sqrt\pi a_\perp}$ 得 $\phi_{0,0}(0)$。
- $n=0$ 散射项 $z\dfrac{e^{ik_zz}}{-ik_z}$：$\partial_z\Big(z\dfrac{e^{ik_zz}}{-ik_z}\Big)\big|_{z\to0}=\dfrac1{-ik_z}=\dfrac{i}{k_z}$，贡献 $\dfrac{\mu g\mathcal R}{\hbar^2\pi a_\perp^2}\dfrac{1}{ik_z}=\dfrac{\gamma\mathcal R}{ik_z}$（注意前面的负号与 $1/(-ik_z)$ 相乘得正）。
- $\mathcal L$ 项：令 $\xi=\dfrac{2z}{a_\perp}$，则 $z\dfrac{a_\perp}{2}\mathcal L=\dfrac{a_\perp^2}{4}\xi\mathcal L(\xi)$，$\partial_z=\dfrac2{a_\perp}\partial_\xi$，故
  $$
  \partial_z\Big[\frac{a_\perp^2}{4}\xi\mathcal L(\xi)\Big]
  =\frac{a_\perp}{2}\partial_\xi\big[\xi\mathcal L(\xi)\big]
  \xrightarrow{\xi\to0}\frac{a_\perp}{2}\left[\frac{2}{\xi}-C+\overline{\mathcal L}\right]_{\rm regular 部分}
  =\frac{a_\perp}{2}\big(-C+\overline{\mathcal L}(\epsilon)\big).
  $$
  这里用 $\xi\mathcal L(\xi)=2+\big(-C+\overline{\mathcal L}\big)\xi+O(\xi^2)$，其导数在 $\xi\to0$ 时趋于常数部分 $-C+\overline{\mathcal L}(\epsilon)$。

汇总（注意 $\mathcal L$ 项前的负号与 $-\dfrac{\mu g\mathcal R}{\hbar^2\pi a_\perp^2}$）：

$$
\mathcal R=\phi_{0,0}(0)+\frac{\gamma\mathcal R}{ik_z}-\frac{\mu g\mathcal R}{\hbar^2\pi a_\perp^2}\frac{a_\perp}{2}\big(-C+\overline{\mathcal L}(\epsilon)\big).
$$

第三项系数 $\dfrac{\mu g\,a_\perp}{2\hbar^2\pi a_\perp^2}=\dfrac{\gamma a_\perp}{2}$，故

$$
\mathcal R=\phi_{0,0}(0)+\frac{\gamma\mathcal R}{ik_z}-\frac{\gamma\mathcal R\,a_\perp}{2}\big(-C+\overline{\mathcal L}(\epsilon)\big).
$$

移项得

$$
\boxed{\;\mathcal R=\frac{\phi_{0,0}(0)}{1+\gamma\left(\dfrac{i}{k_z}+\dfrac{a_\perp}{2}\big[-C+\overline{\mathcal L}(\epsilon)\big]\right)}\;}
$$

---

## 7. 一维散射振幅与一维散射长度

把 $\mathcal R$ 代入 $f_{\rm even}=\dfrac{\gamma\mathcal R}{ik_z\phi_{0,0}(0)}$：

$$
f_{\rm even}=\frac{\gamma}{ik_z}\cdot\frac{1}{1+\gamma\left(\dfrac{i}{k_z}+\dfrac{a_\perp}{2}\big[-C+\overline{\mathcal L}(\epsilon)\big]\right)}
=\frac{\gamma}{ik_z+i\gamma+\dfrac{i\gamma k_z a_\perp}{2}\big[-C+\overline{\mathcal L}(\epsilon)\big]}.
$$

分子分母同乘 $-\dfrac{ik_z}{\gamma}$：

$$
f_{\rm even}=-\frac{1}{1-ik_z\left(\dfrac1\gamma+\dfrac{a_\perp}{2}\big[-C+\overline{\mathcal L}(\epsilon)\big]\right)}.
$$

利用 $\gamma=\dfrac{2a}{a_\perp^2}$，$\dfrac1\gamma=\dfrac{a_\perp^2}{2a}$：

$$
\frac1\gamma+\frac{a_\perp}{2}(-C)=\frac{a_\perp^2}{2a}-\frac{a_\perp C}{2}=-\left[-\frac{a_\perp^2}{2a}+\frac{a_\perp C}{2}\right].
$$

定义一维散射长度

$$
\boxed{\;a_{1D}\equiv-\frac{a_\perp^2}{2a}+\frac{a_\perp C}{2}
=-\frac{a_\perp^2}{2a}\left(1-C\frac{a}{a_\perp}\right)\;}
$$

于是
$$
\frac1\gamma+\frac{a_\perp}{2}\big[-C+\overline{\mathcal L}(\epsilon)\big]=-a_{1D}+\frac{a_\perp}{2}\overline{\mathcal L}(\epsilon).
$$

最终

$$
\boxed{\;f_{\rm even}(k_z)=-\frac{1}{1+ik_z a_{1D}-\dfrac{ik_z a_\perp}{2}\,\overline{\mathcal L}\!\left(-\dfrac{k_z^2a_\perp^2}{4}\right)}\;}
$$

这就是论文的关键结果（式 11）。修正项 $\dfrac{ik_za_\perp}{2}\overline{\mathcal L}=O\big((k_za_\perp)^3\big)$，在低能下可忽略。

> 一维散射长度 $a_{1D}$ 也可用偶宇称相移 $\Delta(k_z)$ 定义为
> $$
> a_{1D}=-\left.\frac{\partial \Delta(k_z)}{\partial k_z}\right|_{k_z\to0^+},
> $$
> 其中偶宇称解 $\Psi\propto\sin\big(k_z|z|+\Delta(k_z)\big)\phi_{0,0}(\rho)$。

---

## 8. 低能极限：有效一维 $\delta$ 势

低能下忽略 $O((k_za_\perp)^3)$ 修正项：

$$
f_{\rm even}(k_z)\approx -\frac{1}{1+ik_z a_{1D}}.
$$

这正是势
$$
U_{1D}(z)=g_{1D}\,\delta(z),\qquad g_{1D}=-\frac{\hbar^2}{\mu a_{1D}},
$$
的偶宇称散射振幅。

把 $a_{1D}$ 代入：

$$
g_{1D}=-\frac{\hbar^2}{\mu a_{1D}}
=\frac{2\hbar^2 a}{\mu a_\perp^2}\left(1-C\frac{a}{a_\perp}\right)^{-1}
=g\,|\phi_{0,0}(0)|^2\left(1-C\frac{a}{a_\perp}\right)^{-1}.
$$

即

$$
\boxed{\;g_{1D}=\frac{2\hbar^2 a}{\mu a_\perp^2}\left(1-C\frac{a}{a_\perp}\right)^{-1}
=\frac{2\pi\hbar^2a}{\mu}\cdot\frac{1}{\pi a_\perp^2}\left(1-C\frac{a}{a_\perp}\right)^{-1}\;}.
$$

其中 $C=1.4603\ldots=-\zeta(1/2)$。

> **适用条件**：紧约束极限 $a_\perp\ll|a|$ 仅在 $s$ 波散射长度 $a$ 远大于真实势的有效程 $r_0$ 时才有意义（$|a|\gg r_0$），因为三维赝势近似仅当 $k\ll r_0^{-1}$ 有效，从而要求 $a_\perp\gtrsim r_0$；二者相容需 $|a|\gg r_0$（近共振情形）。

---

## 9. Tonks 气体与玻色–费米对偶

### 9.1 全反射与不可穿透玻色子

透射系数
$$
T=\big|1+f_{\rm even}+f_{\rm odd}\big|^2=\big|1+f_{\rm even}\big|^2.
$$

低能下 $1+f_{\rm even}=1-\dfrac{1}{1+ik_za_{1D}}=\dfrac{ik_za_{1D}}{1+ik_za_{1D}}$，故

$$
T=\frac{k_z^2a_{1D}^2}{1+k_z^2a_{1D}^2}
\xrightarrow{k_z|a_{1D}|\ll1}0.
$$

即发生**全反射**。这对应一维不可穿透玻色子（Tonks 气体）。形式上它对应无穷强排斥 $g_{1D}\to\infty$。

### 9.2 玻色–费米映射

对周期边界、长 $L$ 的一维盒中 $N$ 个不可穿透玻色子，基态波函数是理想费米气体基态波函数的绝对值：

$$
\Psi_b=\big|\Psi_f\big|,\qquad
\Psi_f(z_1,\ldots,z_N)=\frac{1}{\sqrt{N!\,L^N}}\det\big(e^{ik_j z_{j'}}\big),
$$

$$
k_j=\frac{2\pi}{L}\,j\in\big[-k_{\rm Fermi},\,k_{\rm Fermi}\big],\qquad
k_{\rm Fermi}=\frac{\pi(N-1)}{L}.
$$

其物理来源：一维硬核相互作用使得波函数在粒子坐标重合时为零，与费米子的反对称 Slater 行列式有相同的节点结构，故取绝对值后二者等价。

### 9.3 动量分布与原子数上限

Tonks 气体零温单粒子动量分布 $w(k_z)$ 在原点有峰
$$
w(k_z)\approx\frac{\rho_\infty\sqrt{k_{\rm Fermi}}}{2\pi\sqrt{k_z}},
$$
其中 $\rho_\infty=\pi e^{1/2}2^{-1/3}A^{-6}=0.92418\ldots$，$A=1.2824\ldots$ 为 Glaisher 常数。这是 Tonks 气体的实验标志。

低能条件 $k_{\max}|a_{1D}|\ll1$（$k_{\max}\sim k_{\rm Fermi}$）给出原子数上限

$$
N\ll N^*=\frac{L}{\pi|a_{1D}|}\propto L\,\omega_\perp .
$$

---

## 10. 意义

1. **降维物理的定量基础**：本文第一次严格地（在赝势近似下）把三维冷原子碰撞约化为有效一维 $\delta$ 相互作用，给出了耦合强度 $g_{1D}$ 的解析表达式。后来在原子芯片、光晶格与强拉长势阱中的大量准一维量子气体实验都以此为基础。

2. **约束诱导共振**：$g_{1D}$ 在 $a=a_\perp/C$ 处发散，这就是**约束诱导共振**（CIR），为实验上把一维相互作用从弱到强连续调谐（进而到达 Tonks 极限）提供了机制，成为研究一维多体物理的核心工具。

3. **Tonks 气体的可实现性**：文章预言强横向约束 + 低密度 + 低温条件下可实验实现一维不可穿透玻色子，并给出具体参数（如 $^{87}$Rb、Na）。这一预言在 2004 年由 Bloch、Weiss 等实验组实现，验证了玻色–费米对偶。

4. **玻色–费米对偶与强关联一维体系**：一维硬核玻色子与自由费米子的映射是一维强关联可积系统（Lieb–Liniger 模型、Yang–Gaudin 模型）的极端情形，为理解一维量子液体、自旋–电荷分离等概念提供了实验平台。

---

**一句话总结**：横向简谐约束 + 三维赝势 $\Rightarrow$ 有效一维 $\delta$ 势；其强度 $g_{1D}$ 由三维散射长度 $a$、横向尺寸 $a_\perp$ 与常数 $C=-\zeta(1/2)$ 决定；低能下透射为零，实现 Tonks 气体。
