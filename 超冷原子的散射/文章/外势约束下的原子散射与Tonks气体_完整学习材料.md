# 外势约束中的原子散射与不可穿透玻色子气体：完整学习推导

> **原文**：M. Olshanii, *Atomic Scattering in the Presence of an External Confinement and a Gas of Impenetrable Bosons*, Phys. Rev. Lett. **81**, 938 (1998)；arXiv:cond-mat/9804130。
>
> **本材料的目标与写法**：
>
> 1. 覆盖原文全部要点，但按"物理背景 → 微观模型 → 精确求解 → 低能有效理论 → 多体应用（Tonks 气体）→ 适用范围"的逻辑重新编排。
> 2. 所有公式给出从薛定谔方程或第一性原理出发的推导。**推导顺序严格递进**：后面用到的每个结果都先在前文推导并编号（本材料用 (1)、(2)、… 编号，并标注与原文公式的对应关系）。
> 3. 对理解难点、专业术语与背景知识（冷原子、波导、赝势、散射长度、Green 函数、ζ 函数、Tonks 气体、玻色–费米对偶等）给出专门解释。
> 4. 面向具备量子力学基础（薛定谔方程、谐振子、特殊函数、复变积分）的读者，读完后应能独立复现原文的核心数学计算。

**约定**：两原子质量均为 $m$，约化质量 $\mu=m/2$；取 $\hbar=1$ 的地方会明确说明；二维横向坐标为 $\rho=\sqrt{x^2+y^2}$，纵向坐标为 $z$；横向谐振子频率为 $\omega_\perp$，横向基态尺度为 $a_\perp=\sqrt{\hbar/(\mu\omega_\perp)}$。

---

## 0. 导读

### 0.1 文章做了什么

**问题**：当两个冷原子被一个横向的强简谐势（"波导"或"雪茄形"陷阱）约束、纵向自由运动时，它们在低能下的碰撞有什么规律？能否用一个**有效的一维理论**来替代三维理论？

**答案（原文结果）**：

1. **精确一维散射振幅**（在零程赝势近似下解析求出）：

   $$
   f_{\rm even}(k_z)=-\frac{1}{1+ik_za_{\rm 1D}-\dfrac{ik_za_\perp}{2}\bar{\mathcal L}\!\left(-\dfrac{k_z^2a_\perp^2}{4}\right)},
   \qquad
   a_{\rm 1D}=-\frac{a_\perp^2}{2a}\left(1-\mathcal C\frac{a}{a_\perp}\right),
   $$

   其中 $k_z$ 是纵向（相对）动量，$a$ 是三维 $s$ 波散射长度，$\mathcal C=1.4603\ldots=-\zeta(1/2)$，$\bar{\mathcal L}$ 是一个可解析展开的函数（$=O((k_za_\perp)^2)$，在低能下可忽略）。
2. **有效一维相互作用**：低能下振幅退化为一维 $\delta$ 势的结果

   $$
   f^{\delta}_{\rm even}(k_z)=-\frac{1}{1+ik_za_{\rm 1D}},
   \qquad
   U_{\rm 1D}(z)=g_{\rm 1D}\delta(z),
   \qquad
   g_{\rm 1D}=-\frac{\hbar^2}{\mu a_{\rm 1D}} .
   $$

   横向约束通过 $a_{\rm 1D}$ 和 $g_{\rm 1D}$ 修正了原子间相互作用强度。
3. **约束诱导共振（CIR）**：当 $\mathcal C a/a_\perp=1$（即 $a=a_\perp/\mathcal C\approx0.685\,a_\perp$）时 $g_{\rm 1D}\to\infty$——只需调节三维散射长度（用磁场 Feshbach 共振）就能让一维相互作用从弱变到无穷强。
4. **Tonks 气体**：在低能极限 $k_z|a_{\rm 1D}|\ll1$ 下，一维散射退化为**全反射**，等效于**不可穿透玻色子气体**（Tonks 气体）。该系统与无相互作用的一维费米气体存在一一映射（玻色–费米对偶），基态波函数是费米海 Slater 行列式的绝对值。
5. **实验参数估计**：给出了实现 Tonks 气体的条件与原子数上限 $N\ll N^*=L/(\pi|a_{\rm 1D}|)\propto L\omega_\perp$，以及 $^{87}$Rb、Na 的具体数字（$5\times10^3$ 与 $6\times10^2$）。

**（原文公式编号对照）** 赝势 (1)，薛定谔方程 (2)，$\hat H_\perp$ (3)，能量条件 (4)，渐近形式 (5)，$f_{\rm even},f_{\rm odd}$ (6)，$\eta$ (7)，波函数 (8)，常数 $\mathcal C$ (9)，$\bar{\mathcal L}$ (10)，振幅 (11)，$a_{\rm 1D}$ (12)，$U_{\rm 1D}$ (13)，$g_{\rm 1D}$ (14)，全反射条件 (15)，$\Psi^b=|\Psi^f|$ (16)，$\Psi^f$ (17)，$k_j,k_{\rm Fermi}$ (18)，$N\ll N^*$ (19)。

### 0.2 全篇逻辑链

```
模型：3D 薛定谔方程 + 横向简谐势 + 零程赝势                        §2
   │
   ├─ 两体相对运动分离，约化质量 μ=m/2，a⊥=√(ℏ/μω⊥)              §2.1  (1)(2)(3)
   ├─ 赝势与强度 g=2πℏ²a/μ（Lippmann–Schwinger 推出）             §2.2  (4)–(8)
   ├─ 横向模式：2D 谐振子谱与 |φ_{n,0}(0)|²=1/(πa⊥²)              §2.3  (9)–(13)
   └─ 能量窗口：k_z a⊥ < 2，只有基态道开放                        §2.4  (14)–(17)
   │
按横向模式展开 Ψ=Σ_s ψ_s(z)φ_s(ρ)                                 §3.1  (18)–(26)
   │
   ├─ 相互作用是点源：gδ(r)∂_r(rΨ) → 源强 gφ_s(0)ηδ(z)
   ├─ 1D Green 函数（开放道 outgoing / 闭道 decaying）            §3.2  (27)(28)
   ├─ 求出 ψ_0 与 ψ_s，得到 Ψ(z,0) 的闭合表达式                   §3.3  (29)–(32)
   ├─ 自洽条件 η=∂_z[zΨ]|_{0+} 定出 η，Λ 展开 → C=−ζ(1/2)         §3.4  (33)–(39)
   └─ 精确 1D 振幅 f_even 与 a_1D（含相位定义）                    §3.5  (40)–(43)
   │
低能极限 → δ 势、g_1D、CIR                                        §4    (44)–(49)
   │
强耦合 k_z|a_1D|≪1 → 全反射 → Tonks 气体、玻色–费米对偶          §5    (50)–(57)
   │
适用范围（a⊥≫r0, a⊥≪|a|, |a|≫r0）与总结                         §6    (58)
```

### 0.3 预备知识（如果某一步卡住，先读这里）

| 工具                                                                                                                                           | 用途                          | 本材料位置    |
| ---------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------- | ------------- |
| 两体约化质量与质心分离                                                                                                                         | 从两原子问题得到相对运动方程  | §2.1         |
| 低能散射：散射长度$a$，$s$ 波相移，$f=-1/(1/a+ik)$                                                                                       | 理解赝势为何能替代真实势      | §2.2.1       |
| 三维 Green 函数$-({\mu}/{2\pi\hbar^2})e^{ikr}/r$                                                                                             | 赝势的 Lippmann–Schwinger 解 | §2.2.2       |
| 二维谐振子本征函数                                                                                                                             | 横向模式与$|\phi(0)|^2$     | §2.3、附录 A |
| 一维 Green 函数：$-\frac{i}{2k}e^{ik|z|}$、$-\frac{1}{2\kappa}e^{-\kappa|z|}$                                                              | 解模式方程                    | §3.2、附录 B |
| Riemann/Hurwitz$\zeta$ 函数、Mellin 变换                                        | 提取常数$\mathcal C=-ζ(1/2)$、展开 $\bar{\mathcal L}$ | §3.4、附录 C                 |               |
| 一维$\delta$ 势散射                                                                                                                          | 有效一维理论                  | §4.1         |
| Girardeau 玻色–费米映射                                                                                                                       | Tonks 气体                    | §5.2         |

---

## 1. 物理背景

### 1.1 冷原子、玻色凝聚与原子波导（背景知识）

- 1995 年碱金属原子气体的玻色–爱因斯坦凝聚（BEC）被实现（原文引文 [1]），开启了"超低能量碰撞物理"的新领域：**原子温度极低（nK 量级），德布罗意波长极长，两体碰撞完全由**$s$ 波散射长度 $a$** 决定。**
- 同时，人们可以用**光诱导**[2]或**磁场诱导**[3]的"原子波导（atom waveguide）"把原子限制在准一维的几何中。**波导**在这里指：一个方向（$z$）让原子自由运动，另外两个方向（$x,y$）用势场强烈约束，形成一个细长的"雪茄"。
- 当横向约束很强、温度足够低时，原子的横向运动被冻结在谐振子基态，物理变为**有效一维（quasi-1D）**的。一维气体与三维气体有本质区别：一维中两粒子无法"绕开"彼此，相互作用效应被大大增强，而且可以出现只存在于一维的奇特现象（如自旋–电荷分离、玻色–费米对偶等）。

**术语**：*准一维*指体系的横向尺寸远小于纵向尺度、且横向激发被能量禁戒；*雪茄形（cigar-shaped）陷阱*指各向异性很大的谐振子陷阱 $\omega_\perp\gg\omega_z$ 的极限（本文进一步取 $\omega_z\to0$）。

### 1.2 核心问题

原子间的"真实"相互作用是三维的短程势（作用范围 $\sim$ 几十个玻尔半径），由 $s$ 波散射长度 $a$ 描述。问题是：

1. 在横向约束下，两个原子的**有效一维碰撞**由什么参数描述？
2. 能否把三维问题严格约化为一个**一维 $\delta$ 势**问题？其耦合常数是多少？
3. 低能极限下会发生什么？能否人为调节相互作用强度？
4. 这些结果对**多体**（N 原子）体系意味着什么？

本文在**零程赝势（pseudopotential）近似**下解析地回答了这些问题。赝势近似的含义与条件见 §2.2 与 §6.1。

### 1.3 文章的结论预览（后文会逐一推导）

- 约束并不只是"把三维问题降维"：横向的**虚激发**（碰撞过程中暂时跃迁到高横向模式）会对低能散射产生一个**修正**，其后果是把 $a$ 换成 $a_{\rm 1D}=-\frac{a_\perp^2}{2a}(1-\mathcal C a/a_\perp)$。这个修正可以在实验上被"共振"放大（CIR）。
- 低能极限（$k_z|a_{\rm 1D}|\ll1$）下，一维透射系数趋于零——系统的两体散射等价于**硬核**（不可穿透）。因此稀薄的准一维玻色气体成为 **Tonks 气体**：一个强关联的玻色系统，其性质与理想费米气体一一对应（Girardeau 1960）。
- Tonks 气体的动量分布在 $k_z\to0$ 出现 $1/\sqrt{k_z}$ 发散，可作为实验判据。

---

## 2. 微观模型：从三维薛定谔方程出发

### 2.1 两原子、质心与相对运动（推导）

考虑两个质量均为 $m$ 的原子，三维坐标 $\mathbf r_1,\mathbf r_2$。它们在横向简谐陷阱中，陷阱频率为 $\omega_\perp$；相互作用是短程两体势 $V(\mathbf r_1-\mathbf r_2)$。总哈密顿量（相对运动部分我们最终只关心 $z$ 与 $\rho$）：

$$
H=\frac{\mathbf p_1^2}{2m}+\frac{\mathbf p_2^2}{2m}
+\frac12 m\omega_\perp^2\big(\rho_1^2+\rho_2^2\big)+V(\mathbf r_1-\mathbf r_2)
+\underbrace{\Big[\frac12 m\omega_z^2(z_1^2+z_2^2)\Big]}_{\text{纵向约束（本文取 }\omega_z\to0\text{）}},
\tag{1}
$$

其中 $\rho_i^2=x_i^2+y_i^2$。引入质心与相对坐标

$$
\mathbf R=\frac{\mathbf r_1+\mathbf r_2}{2},\qquad
\mathbf r=\mathbf r_2-\mathbf r_1,
\qquad
\mathbf P=\mathbf p_1+\mathbf p_2,\qquad
\mathbf p=\frac{\mathbf p_2-\mathbf p_1}{2},
$$

**推导（动能与横向势的分离）**：由 $\mathbf r_1=\mathbf R-\mathbf r/2$、$\mathbf r_2=\mathbf R+\mathbf r/2$ 得

$$
\mathbf p_1=\frac{\mathbf P}{2}-\mathbf p,\qquad
\mathbf p_2=\frac{\mathbf P}{2}+\mathbf p,
$$

于是

$$
\frac{\mathbf p_1^2+\mathbf p_2^2}{2m}
=\frac{\mathbf P^2}{2(2m)}+\frac{\mathbf p^2}{2(m/2)}
=\frac{\mathbf P^2}{2M}+\frac{\mathbf p^2}{2\mu},
\qquad M=2m,\ \ \mu=\frac m2 .
$$

对横向势：

$$
\rho_1^2+\rho_2^2=2\rho_R^2+\frac{\rho^2}{2},
$$

所以

$$
\frac12 m\omega_\perp^2(\rho_1^2+\rho_2^2)
=\frac12 M\omega_\perp^2\rho_R^2+\frac12\mu\omega_\perp^2\rho^2 .
$$

质心运动与相对运动完全分离。我们只研究相对运动，其哈密顿量为

$$
\boxed{\;
H_{\rm rel}=\frac{\hat p_z^2+\hat p_x^2+\hat p_y^2}{2\mu}
+\frac12\mu\omega_\perp^2(x^2+y^2)
+V(\mathbf r)\;}
\tag{2}
$$

这正是原文式 (2)(3) 的相对运动部分（原文把 $\hat p_z^2/2\mu$ 与横向部分分开写）。**要点**：横向陷阱频率对相对运动不变（还是 $\omega_\perp$），但振子质量是约化质量 $\mu=m/2$，因此相对运动的横向振子长度

$$
\boxed{\;a_\perp\equiv\sqrt{\frac{\hbar}{\mu\omega_\perp}}=\sqrt{\frac{2\hbar}{m\omega_\perp}}\;}
\tag{3}
$$

才是后文反复出现的"横向尺寸"。$a_\perp$ 也是横向基态波函数（高斯）的宽度（见 §2.3）。

### 2.2 短程相互作用的低能描述：Huang 赝势

真实原子间势 $V(r)$ 的力程很短（记 $r_0$ 为有效力程），在本文关心的低能下，只需一个参数 $s$ 波散射长度 $a$ 就可描述两体散射。为了能把问题解析求解，原文把 $V$ 替换为**零程赝势**（原文式 (1)）：

$$
\boxed{\;
U(\mathbf r)=g\,\delta(\mathbf r)\frac{\partial}{\partial r}\big(r\,\cdot\,\big)\;}
\tag{4}
$$

本节推导它的来历、正规化算符的必要性以及耦合常数

$$
\boxed{\;g=\frac{2\pi\hbar^2 a}{\mu}\;}
\tag{5}
$$

（原文式 (1)；注意原文用 $\mu=m/2$，所以等价于常见的 $4\pi\hbar^2a/m$。）

#### 2.2.1 低能散射与散射长度（背景推导）

**三维 $s$ 波散射**：把相对运动的散射态按分波展开，低能下只有 $L=0$ 重要，其径向波函数 $u(r)=r\psi(r)$ 满足

$$
-\frac{\hbar^2}{2\mu}u''+V(r)u=Eu,\qquad E=\frac{\hbar^2 k^2}{2\mu},
$$

在 $r>r_0$（势外）$u(r)\propto\sin(kr+\delta_0)$。**散射长度**的定义来自零能动量极限的相移：

$$
k\cot\delta_0(k)\ \xrightarrow{k\to0}\ -\frac1a
\qquad\Longleftrightarrow\qquad
\delta_0(k\to0)\simeq -ka .
$$

分波散射振幅为

$$
f=\frac{e^{2i\delta_0}-1}{2ik}=\frac{1}{k\cot\delta_0-ik}
\ \xrightarrow{\text{低能}}\ -\frac{1}{1/a+ik}
\ \xrightarrow{k\to0}\ -a .
\tag{6}
$$

**推导（第二步）**：$\dfrac{e^{2i\delta}-1}{2ik}=\dfrac{e^{i\delta}\sin\delta}{k}=\dfrac{1}{k\cot\delta-ik}$。∎

**物理意义**：$a$ 是零能量散射振幅的相反数，且零能波函数在势外为 $u(r)\propto r-a$，即 $\psi\propto(1-a/r)$。

#### 2.2.2 为什么是 $\delta$ 函数：赝势的形式（背景）

零程势的思想：既然低能散射只记得一个参数 $a$，就构造一个作用范围为**零**、但强度与 $a$ 相关的势，使它给出相同的低能散射振幅。裸 $\delta$ 势 $g\delta(\mathbf r)$ 的散射解在原点附近行为不对（含有未正规化的 $1/r$），因此 Huang 引入

$$
U=g\,\delta(\mathbf r)\frac{\partial}{\partial r}(r\,\cdot\,),
$$

其中算符 $\frac{\partial}{\partial r}(r\,\cdot\,)$ 先取 $r\psi$ 并求导，再取 $r\to0$ 的**正则部分**；它自动"减掉" $1/r$ 发散项。下面用 Lippmann–Schwinger 方程从这个势**推出** $g$ 与 $a$ 的关系 (5)。

#### 2.2.3 从薛定谔方程推出 $g=2\pi\hbar^2a/\mu$

在自由空间（没有横向势）中考虑势 (4) 的散射。三维自由 Green 函数（出射波边界条件）为

$$
G_0(\mathbf r,\mathbf r')=-\frac{\mu}{2\pi\hbar^2}\frac{e^{ik|\mathbf r-\mathbf r'|}}{|\mathbf r-\mathbf r'|},
\tag{7}
$$

**验证**：$(∇^2+k^2)\frac{e^{ikr}}{r}=-4\pi\delta^{(3)}(\mathbf r)$，故 $G_0$ 满足 $(\hat H_0-E)G_0=\delta^{(3)}(\mathbf r-\mathbf r')$。∎

**推导**：Lippmann–Schwinger 方程为

$$
\psi(\mathbf r)=e^{ikz}+\int d^3r'\,G_0(\mathbf r,\mathbf r')\,U(\mathbf r')\psi(\mathbf r').
$$

注意到势 (4) 中的 $\delta(\mathbf r')$ 把积分变成在原点取值，并定义

$$
\eta\equiv\left.\frac{\partial}{\partial r'}\big(r'\psi\big)\right|_{r'\to0},
$$

则

$$
\int d^3r'\,G_0(\mathbf r,\mathbf r')U(\mathbf r')\psi(\mathbf r')
=g\,\eta\int d^3r'\,G_0(\mathbf r,\mathbf r')\delta^{(3)}(\mathbf r')
=-\,g\eta\,\frac{\mu}{2\pi\hbar^2}\frac{e^{ikr}}{r}.
$$

于是

$$
\psi(\mathbf r)=e^{ikz}-\frac{\mu g}{2\pi\hbar^2}\,\eta\,\frac{e^{ikr}}{r}
=e^{ikz}+f\,\frac{e^{ikr}}{r},
\qquad
f=-\frac{\mu g}{2\pi\hbar^2}\eta .
\tag{8}
$$

**自洽条件**：$\eta$ 的定义要求它是波函数在原点附近的**正则部分**。把 $\psi$ 在 $r\to0$ 展开：

$$
e^{ikz}=1+ikz+\cdots,\qquad
\frac{e^{ikr}}{r}=\frac1r+ik+O(r),
$$

所以

$$
r\psi=\big(1+ikz+\cdots\big)r+f\,e^{ikr}
\ \Longrightarrow\
\eta=\left.\frac{\partial (r\psi)}{\partial r}\right|_{r\to0}
=1+ikf+O(r).
$$

再利用 (8) 得

$$
\eta=1+ik\Big(-\frac{\mu g}{2\pi\hbar^2}\eta\Big)
\ \Longrightarrow\
\eta=\frac{1}{1+i k\,\frac{\mu g}{2\pi\hbar^2}} .
$$

代入 (8)：

$$
f=-\frac{\mu g}{2\pi\hbar^2}\cdot\frac{1}{1+i k\,\frac{\mu g}{2\pi\hbar^2}}
=-\frac{1}{\frac{2\pi\hbar^2}{\mu g}+ik}.
$$

与低能标准结果 (6) $f=-1/(1/a+ik)$ 比较，得到

$$
\boxed{\;\frac{2\pi\hbar^2}{\mu g}=\frac1a\ \Longleftrightarrow\ g=\frac{2\pi\hbar^2a}{\mu}\;}
\qquad\text{(即式 (5))}
$$

**结论**：在低能、$s$ 波条件下，赝势 (4) 与真实势给出相同的散射长度 $a$。赝势近似的**适用条件**是 $k\ll r_0^{-1}$（见 §6.1）。∎

### 2.3 横向束缚：二维谐振子（推导）

相对运动的横向哈密顿量为

$$
\hat H_\perp=\frac{\hat p_x^2+\hat p_y^2}{2\mu}+\frac12\mu\omega_\perp^2(x^2+y^2).
\tag{9}
$$

**本征值与基态**（推导见附录 A）：在极坐标 $(\rho,\varphi)$ 中，

$$
\hat H_\perp\phi_{n,m}(\rho,\varphi)=E_n\phi_{n,m},\qquad
\boxed{\;E_n=\hbar\omega_\perp(n+1),\qquad n=0,1,2,\ldots\;}
\tag{10}
$$

（二维各向同性谐振子的能谱是 $n_x+n_y+1$ 的简并谱）。$m$ 是角动量量子数，$m=0,\pm2,\ldots$（$n$ 偶）或 $m=\pm1,\pm3,\ldots$（$n$ 奇）。基态为

$$
\boxed{\;\phi_{0,0}(\rho)=\frac{1}{\sqrt{\pi}\,a_\perp}e^{-\rho^2/(2a_\perp^2)},\qquad
a_\perp=\sqrt{\frac{\hbar}{\mu\omega_\perp}}\;}
\tag{11}
$$

归一化：$\int|\phi_{0,0}|^2d^2\rho=\frac{1}{\pi a_\perp^2}\int_0^\infty e^{-\rho^2/a_\perp^2}2\pi\rho\,d\rho=1$。因此在 $\rho=0$ 处

$$
\phi_{0,0}(0)=\frac{1}{\sqrt{\pi}a_\perp},\qquad
|\phi_{0,0}(0)|^2=\frac{1}{\pi a_\perp^2}.
\tag{12}
$$

**关键性质（后文反复使用）**：对所有角动量为零的态（$m=0$，即 $n=0,2,4,\ldots$），

$$
\boxed{\;|\phi_{n,0}(0)|^2=\frac{1}{\pi a_\perp^2}\qquad(n=0,2,4,\ldots)\;}
\tag{13}
$$

**推导**：$m=0$ 的径向波函数为（详细推导与解的性质见附录 A）
$\phi_{n,0}(\rho)=\frac{1}{\sqrt{2\pi}}R_{n,0}(\rho)$，
$R_{n,0}(\rho)=\sqrt{\frac{2}{a_\perp^2}}\,e^{-\rho^2/(2a_\perp^2)}\,L_{n/2}\!\big(\frac{\rho^2}{a_\perp^2}\big)$，
其中 $L_k$ 是拉盖尔多项式、$L_k(0)=1$（附录 A）。于是在原点

$$
\phi_{n,0}(0)=\frac{1}{\sqrt{2\pi}}\sqrt{\frac{2}{a_\perp^2}}=\frac{1}{\sqrt{\pi}a_\perp}\qquad(\forall n\ \text{偶}),
$$

与 $n$ 无关。∎

### 2.4 能量窗口与散射渐近形式

#### 2.4.1 只有基态横向模式开放（推导）

考虑相对运动的散射态，总能量为

$$
E=\frac{\hbar^2k_z^2}{2\mu}+\hbar\omega_\perp
\qquad(k_z:\text{纵向相对动量}),
\tag{14}
$$

即入射波处于横向基态 $\phi_{0,0}$（原文假设 (i)）。

**能量约束**：若 $E$ 超过某个激发态的阈值，碰撞就可以把横向模式激发，散射变成非弹性的，一维图像失效。由角动量守恒，轴对称为零的相互作用只能耦合到 $m=0$ 的态；而 $m=0$ 的最低激发态是 $n=2$（$E_{2,0}=3\hbar\omega_\perp$）。因此要求

$$
\frac{\hbar^2k_z^2}{2\mu}<E_{n=2,m=0}-E_{n=0,m=0}=2\hbar\omega_\perp .
\tag{15}
$$

用 (3) 改写：

$$
\boxed{\;k_z^2a_\perp^2<4\;}
\tag{16}
$$

（原文式 (4)）。注意：(i) 更高 $n$ 的 $m=0$ 态在**碰撞过程中**作为虚中间态仍会被激发（这正是下节模式展开要处理的核心）；(ii) 条件 (16) 保证所有 $n\ge2$ 的模式都是"闭合道"（能量上不可及），它们的波函数在 $|z|\to\infty$ 处指数衰减。

#### 2.4.2 渐近形式与一维散射振幅（推导）

在 $|z|\to\infty$，横向只能处于基态 $\phi_{0,0}$（激发模式的波函数指数衰减已消失）。把散射态分解为关于 $z\to-z$ 的偶、奇部分。入射波 $e^{ik_zz}\phi_{0,0}$ 可拆成 $\cos(k_zz)\phi_{0,0}$（偶）与 $i\sin(k_zz)\phi_{0,0}$（奇）；相互作用是 $z$ 的偶函数（关于 $z=0$ 对称），因此偶、奇道互不耦合。

- **奇道**：自由解 $i\sin(k_zz)$ 在 $z=0$ 处为零（远离相互作用点），所以奇道完全不受影响，$f_{\rm odd}=0$。
- **偶道**：散射解为 $\cos(k_zz)+f_{\rm even}e^{ik_z|z|}$；乘上平面波分量后，总的渐近形式为

$$
\boxed{\;
\Psi(z,\rho)\ \xrightarrow{|z|\to\infty}\
\Big\{e^{ik_zz}+f_{\rm even}(k_z)\,e^{ik_z|z|}
+f_{\rm odd}(k_z)\,{\rm sign}(z)\,e^{ik_z|z|}\Big\}\phi_{0,0}(\rho)\;}
\tag{17}
$$

（原文式 (5)），且已论证 $f_{\rm odd}=0$。$f_{\rm even}$ 称为**一维散射振幅**（偶宇称）。本文的主要任务就是求出 $f_{\rm even}(k_z)$。**注意**：式 (17) 用 $|z|$ 与 $\mathrm{sign}(z)$ 表示"出射波"；在一维中，出射波是 $e^{\pm ik_z|z|}$。

---

## 3. 精确求解：横向模式展开 + Green 函数

本节是全文的技术核心：把三维方程投影到横向本征态上，得到无穷多个一维耦合方程；利用"相互作用是原点处的点源"这一事实，用 Green 函数逐道求解，最后用一个自洽条件定出振幅。

### 3.1 横向模式展开与投影（推导）

把散射态按横向本征态展开。由 §2.4 的讨论，只需保留 $m=0$ 的态。记 $n=2s\ (s=0,1,2,\ldots)$，并简写

$$
\phi_s(\rho)\equiv\phi_{2s,0}(\rho),\qquad
E_s=\hbar\omega_\perp(2s+1),\qquad
\phi_s(0)=\frac{1}{\sqrt{\pi}a_\perp}\ \ (\text{由式 (13)}).
\tag{18}
$$

展开式：

$$
\boxed{\;\Psi(z,\rho)=\sum_{s=0}^{\infty}\psi_s(z)\,\phi_s(\rho)\;}
\tag{19}
$$

（原文在 (8) 式中只写出 $\rho=0$ 的形式，但推导基于此展开。）

**推导（投影方程）**：把 (19) 代入相对运动薛定谔方程

$$
\left[\frac{\hat p_z^2}{2\mu}+g\,\delta(\mathbf r)\frac{\partial}{\partial r}(r\,\cdot\,)+\hat H_\perp\right]\Psi=E\Psi ,
\tag{20}
$$

（原文式 (2)），并注意

$$
\delta(\mathbf r)\,\frac{\partial}{\partial r}(r\Psi)
=\delta(\rho)\delta(z)\,\eta\equiv\delta(\mathbf r)\,\eta,
\qquad
\eta\equiv\left.\frac{\partial}{\partial r}\big(r\Psi\big)\right|_{r\to0},
\tag{21}
$$

这里用了 $\delta(z)h(z)=h(0)\delta(z)$、$h(0)=\eta$（$h(z)$ 在 $z=0$ 正则），这正是赝势正规化算符的关键性质（原文式 (7)：

$$
\eta=\left.\frac{\partial}{\partial r}\big[r\Psi(\mathbf r)\big]\right|_{r\to0}
=\left.\frac{\partial}{\partial z}\big[z\Psi(z,\rho=0)\big]\right|_{z\to0+}.
\tag{22}
$$

）。用 $\int d^2\rho\,\phi_s^*(\rho)$ 投影 (20)：

$$
-\frac{\hbar^2}{2\mu}\psi_s''(z)+E_s\psi_s(z)+g\,\phi^*_s(0)\,\eta\,\delta(z)=E\,\psi_s(z),
$$

即

$$
\boxed{\;\psi_s''(z)-\kappa_s^2\,\psi_s(z)
=+\frac{2\mu g}{\hbar^2}\,\phi_s^*(0)\,\eta\,\delta(z)\;}
\tag{23}
$$

其中

$$
\kappa_s^2\equiv\frac{2\mu}{\hbar^2}\big(E_s-E\big).
\tag{24}
$$

**推导（$\kappa_s$ 的显式形式）**：把 (14)、(18) 代入 (24)：

- $s=0$（开放道）：$E_0-E=\hbar\omega_\perp-\hbar\omega_\perp-\frac{\hbar^2k_z^2}{2\mu}=-\frac{\hbar^2k_z^2}{2\mu}$，故

$$
\boxed{\;\kappa_0=ik_z\;}
\tag{25}
$$

- $s\ge1$（闭合道）：$E_s-E=\hbar\omega_\perp(2s+1)-\hbar\omega_\perp-\frac{\hbar^2k_z^2}{2\mu}=2s\hbar\omega_\perp-\frac{\hbar^2k_z^2}{2\mu}$，故

$$
\kappa_s^2=\frac{2\mu}{\hbar^2}\Big(2s\hbar\omega_\perp-\frac{\hbar^2k_z^2}{2\mu}\Big)
=\frac{4\mu\omega_\perp}{\hbar}s-k_z^2
=\frac{4}{a_\perp^2}\Big(s+\epsilon\Big),
\qquad
\boxed{\;\epsilon\equiv-\frac{k_z^2a_\perp^2}{4}\;}
\tag{26}
$$

（用了 (3)：$4\mu\omega_\perp/\hbar=4/a_\perp^2$）。由 (16) 有 $\epsilon>-1$，从而 $\kappa_s>0$（$s\ge1$）：闭合道波函数随 $|z|$ 指数衰减，这正是"虚激发"的数学体现。∎

**难点的澄清**：方程 (23) 中右端是 $\delta(z)$ 点源，强度 $\propto\phi_s^*(0)\eta$。这不是额外假设，而是赝势 $\delta(\mathbf r)$ 的严格结果：$\delta(\rho)\delta(z)\eta$ 投影后就是 $\phi_s^*(0)\eta\delta(z)$。$\eta$ 是唯一的未知常数，它承担了"原点处相互作用有多强"的全部信息，由 (22) 自洽确定。

### 3.2 一维 Green 函数（推导）

解方程 (23) 需要两个一维 Green 函数。

**(a) 开放道（出射波）**：$G_0(z)$ 满足

$$
\big(\partial_z^2+k_z^2\big)G_0(z)=\delta(z),
\qquad
G_0(z)\ \text{在}\ |z|\to\infty\ \text{为出射波}.
$$

**推导**：作 Fourier 变换 $G_0(z)=\int\frac{dq}{2\pi}\tilde G_0(q)e^{iqz}$，得 $\tilde G_0(q)=1/(k_z^2-q^2+i0^+)$，极点 $q=\pm(k_z+i0^+)$。对 $z>0$ 在上半平面闭合围道，拾取 $q=+k_z$ 的极点；对 $z<0$ 在下半平面拾取 $q=-k_z$。结果是

$$
\boxed{\;G_0(z)=-\frac{i}{2k_z}e^{ik_z|z|}\;}
\tag{27}
$$

**验证跳变条件**：$G_0'(0^+)-G_0'(0^-)=1$（对 $(\partial_z^2+k_z^2)$ 积分 $\int_{0^-}^{0^+}dz$ 得 $G_0'(0^+)-G_0'(0^-)=\int\delta=1$）。由 (27)：$G_0'=-\frac{i}{2k_z}\cdot ik_z\,\mathrm{sign}(z)e^{ik_z|z|}=\frac12\mathrm{sign}(z)e^{ik_z|z|}$，故跳变为 $1$ ✓。∎

**(b) 闭合道（指数衰减）**：$G_s(z)$ 满足 $(\partial_z^2-\kappa_s^2)G_s=\delta(z)$，且 $|z|\to\infty$ 衰减：

$$
\boxed{\;G_s(z)=-\frac{1}{2\kappa_s}e^{-\kappa_s|z|}\;}
\tag{28}
$$

**验证**：同样由跳变条件：$G_s'=\frac12\mathrm{sign}(z)e^{-\kappa_s|z|}$，跳变 $1$ ✓；且 $(\partial_z^2-\kappa_s^2)G=0$（$z\neq0$）。∎

### 3.3 各通道波函数与 $\Psi(z,0)$

#### 3.3.1 开放道

由 (23)(25)，$s=0$ 的方程与入射条件为

$$
\big(\partial_z^2+k_z^2\big)\psi_0
=\frac{2\mu g}{\hbar^2}\phi_0^*(0)\eta\,\delta(z),
\qquad
\psi_0\ \xrightarrow{z\to-\infty}\ e^{ik_zz}.
$$

用 (27)：

$$
\boxed{\;\psi_0(z)=e^{ik_zz}
+\frac{2\mu g}{\hbar^2}\phi_0^*(0)\eta\,G_0(z)
=e^{ik_zz}
-\frac{i\mu g}{\hbar^2k_z}\phi_0^*(0)\eta\,e^{ik_z|z|}\;}
\tag{29}
$$

#### 3.3.2 闭合道

由 (23)(26)，$s\ge1$：

$$
\big(\partial_z^2-\kappa_s^2\big)\psi_s
=\frac{2\mu g}{\hbar^2}\phi_s^*(0)\eta\,\delta(z),
\qquad
\psi_s\ \xrightarrow{|z|\to\infty}\ 0 .
$$

用 (28)：

$$
\boxed{\;\psi_s(z)
=-\frac{\mu g}{\hbar^2}\frac{\phi_s^*(0)\eta}{\kappa_s}\,e^{-\kappa_s|z|},
\qquad s\ge1\;}
\tag{30}
$$

#### 3.3.3 在 $\rho=0$ 处的总波函数（推导）

把 (18)(25)(26)(29)(30) 代入 (19)，取 $\rho=0$（只有 $\phi_s(0)$ 出现）：

$$
\Psi(z,0)=\sum_{s=0}^\infty\psi_s(z)\phi_s(0)
=\underbrace{\phi_0(0)e^{ik_zz}}_{(\rm I)}
\underbrace{-\frac{i\mu g}{\hbar^2k_z}\big|\phi_0(0)\big|^2\eta\,e^{ik_z|z|}}_{(\rm II)}
\underbrace{-\frac{\mu g}{\hbar^2}\sum_{s=1}^\infty\frac{|\phi_s(0)|^2\eta}{\kappa_s}e^{-\kappa_s|z|}}_{(\rm III)} .
$$

利用 (12)(13)：$|\phi_s(0)|^2=1/(\pi a_\perp^2)$；并用 (3) 与 (5) 化简系数。直接代入 $g=2\pi\hbar^2a/\mu$：

$$
\frac{\mu g}{\hbar^2}=\frac{\mu}{\hbar^2}\cdot\frac{2\pi\hbar^2a}{\mu}=2\pi a,
\qquad
\frac{\mu g}{\hbar^2 k_z}\big|\phi_0(0)\big|^2
=\frac{2\pi a}{k_z}\cdot\frac{1}{\pi a_\perp^2}
=\frac{2a}{k_z a_\perp^2},
$$

以及 $\kappa_s=\frac{2}{a_\perp}\sqrt{s+\epsilon}$。于是

$$
\boxed{\;
\Psi(z,0)=
\frac{1}{\sqrt{\pi}a_\perp}e^{ik_zz}
-\frac{ig\mu\eta}{\pi\hbar^2k_za_\perp^2}e^{ik_z|z|}
-\frac{g\mu\eta}{2\pi\hbar^2a_\perp}\,
\Lambda\!\left[\frac{2|z|}{a_\perp},\,\epsilon\right]
\;}
\tag{31}
$$

（原文式 (8)），其中

$$
\boxed{\;
\Lambda[\xi,\epsilon]\equiv\sum_{s'=1}^{\infty}
\frac{e^{-\sqrt{s'+\epsilon}\,\xi}}{\sqrt{s'+\epsilon}}
\;}
\tag{32}
$$

（原文在求和指标上写 $s'=n/2$，本文统一写成 $s'=1,2,3,\ldots$，即 $n=2,4,6,\ldots$）。

**解读**：第 (I) 项是入射波；第 (II) 项是开放道的出射（散射）波；第 (III) 项是碰撞过程中**虚激发的横向云**——它对 $|z|\gtrsim a_\perp$ 指数衰减，但在碰撞区（$|z|\lesssim a_\perp$）有不可忽略的贡献，正是它产生了对 $a$ 的修正。

### 3.4 自洽条件与 $\Lambda$ 函数的展开

#### 3.4.1 用 $\eta$ 的定义（推导）

把 (31) 代入 (22)：

$$
\eta=\left.\frac{\partial}{\partial z}\big[z\,\Psi(z,0)\big]\right|_{z\to0+}.
$$

对 (31) 的三项分别计算上式右端的导数。第 (I) 项：

$$
z\cdot\frac{1}{\sqrt{\pi}a_\perp}e^{ik_zz}
\ \Longrightarrow\
\frac{\partial}{\partial z}\Big|_{0^+}=\frac{1}{\sqrt{\pi}a_\perp}.
$$

对第 (II) 项（在 $z>0$ 处 $e^{ik_z|z|}=e^{ik_zz}$）：

$$
z\cdot\Big(-\frac{ig\mu\eta}{\pi\hbar^2k_za_\perp^2}e^{ik_zz}\Big)
\ \Longrightarrow\
\frac{\partial}{\partial z}\Big|_{0^+}=-\frac{ig\mu\eta}{\pi\hbar^2k_za_\perp^2}.
$$

对第 (III) 项：需要 $\Lambda$ 在 $\xi\to0$ 的展开。**关键点**（原文特别提醒）：不能把 $\partial_z$ 与求和交换，因为 $\Lambda$ 在 $\xi\to0$ 时不一致收敛（$\Lambda\sim2/\xi$ 发散）。用后文 (35) 的展开

$$
\Lambda[\xi,\epsilon]=\frac{2}{\xi}+\mathcal L(\epsilon)+\mathcal L_1(\epsilon)\xi+\ldots,
\qquad
\mathcal L(\epsilon)=-\mathcal C+\bar{\mathcal L}(\epsilon),
\tag{33}
$$

则在 $z>0$、$\xi=2z/a_\perp\to0$ 时

$$
z\,\Lambda\!\left[\frac{2z}{a_\perp},\epsilon\right]
=z\Big[\frac{a_\perp}{z}+\mathcal L(\epsilon)+O(z)\Big]
=a_\perp+z\,\mathcal L(\epsilon)+O(z^2),
$$

所以

$$
\frac{\partial}{\partial z}\Big[z\Lambda\Big]\Big|_{0^+}=\mathcal L(\epsilon)=-\mathcal C+\bar{\mathcal L}(\epsilon).
$$

于是第 (III) 项的贡献为 $-\dfrac{g\mu\eta}{2\pi\hbar^2a_\perp}\big(-\mathcal C+\bar{\mathcal L}\big)$。把三项相加：

$$
\eta=\frac{1}{\sqrt{\pi}a_\perp}
-\frac{ig\mu\eta}{\pi\hbar^2k_za_\perp^2}
-\frac{g\mu\eta}{2\pi\hbar^2a_\perp}\big(-\mathcal C+\bar{\mathcal L}\big).
\tag{34}
$$

#### 3.4.2 $\Lambda$ 展开、常数 $\mathcal C$ 与函数 $\bar{\mathcal L}$（推导）

**推导（Euler–Maclaurin/解析延拓法）**：在 $\Lambda$ 中"减掉再加"积分表示

$$
\sum_{s'=1}^\infty\int_{s'-1}^{s'}\frac{ds''}{\sqrt{s''}}e^{-\sqrt{s''}\xi}
=\int_0^\infty\frac{ds''}{\sqrt{s''}}e^{-\sqrt{s''}\xi}
=\frac{2}{\xi} .
\tag{35}
$$

于是

$$
\Lambda[\xi,\epsilon]-\frac{2}{\xi}
=\sum_{s'=1}^{\infty}\left[
\frac{e^{-\sqrt{s'+\epsilon}\xi}}{\sqrt{s'+\epsilon}}
-\int_{s'-1}^{s'}ds''\frac{e^{-\sqrt{s''}\xi}}{\sqrt{s''}}
\right].
$$

该级数在 $\xi\to0$ 时一致收敛（原文注），因此可以令 $\xi=0$：

$$
\mathcal L(\epsilon)\equiv\lim_{\xi\to0}\left(\Lambda-\frac2\xi\right)
=\sum_{s'=1}^{\infty}\left[
\frac{1}{\sqrt{s'+\epsilon}}
-2\big(\sqrt{s'}-\sqrt{s'-1}\big)
\right].
\tag{36}
$$

令 $\epsilon=0$ 并用 Riemann $\zeta$ 函数的解析延拓可以证明（附录 C；原文式 (9) 的另一种等价写法）

$$
\boxed{\;
\mathcal L(0)=\sum_{s'=1}^{\infty}\left[
\frac{1}{\sqrt{s'}}-2\big(\sqrt{s'}-\sqrt{s'-1}\big)
\right]
=\zeta\!\Big(\frac12\Big)=-\mathcal C,
\qquad
\mathcal C=1.4603\ldots
\;}
\tag{37}
$$

（更标准的写法：$\mathcal C=\lim_{s\to\infty}\big[\int_0^s\frac{ds'}{\sqrt{s'}}-\sum_{s'=1}^s\frac{1}{\sqrt{s'}}\big]$，原文式 (9)，且 $\mathcal C=-\zeta(1/2)$。）

**对 $\epsilon$ 展开（推导 $\bar{\mathcal L}$）**：把 (36) 中的 $1/\sqrt{s'+\epsilon}$ 在 $\epsilon=0$ 处展开：

$$
\frac{1}{\sqrt{s'+\epsilon}}
=\frac{1}{\sqrt{s'}}\sum_{n=0}^\infty
(-1)^n\frac{(2n-1)!!}{2^nn!}\Big(\frac{\epsilon}{s'}\Big)^n
=\frac{1}{\sqrt{s'}}
+\sum_{n=1}^\infty(-1)^n\frac{(2n-1)!!}{2^nn!}\frac{\epsilon^n}{s'^{\,n+1/2}} .
$$

代入 (36)，利用 $\sum_{s'=1}^\infty s'^{-(n+1/2)}=\zeta(n+1/2)$，并减掉 $\epsilon=0$ 的常数 (37)：

$$
\boxed{\;
\mathcal L(\epsilon)=-\mathcal C+\bar{\mathcal L}(\epsilon),
\qquad
\bar{\mathcal L}(\epsilon)=\sum_{n=1}^{\infty}(-1)^n\,
\frac{\zeta\!\big(\frac{1+2n}{2}\big)(2n-1)!!}{2^nn!}\,\epsilon^n
\;}
\tag{38}
$$

（原文式 (10)）。∎

**几点说明**：

- $\bar{\mathcal L}(0)=0$，故 $\mathcal L(0)=-\mathcal C$。
- 第一项 $n=1$：$-{\zeta(3/2)}\,\epsilon/2$，即 $\bar{\mathcal L}=O(k_z^2a_\perp^2)$。
- 附录 C 给出**更优雅的等价推导**：用 Mellin 变换可得 $\mathcal L(\epsilon)=\zeta(1/2,1+\epsilon)$（Hurwitz $\zeta$ 函数），它直接给出 (37)(38) 两式。

### 3.5 一维散射振幅与一维散射长度（推导）

**第一步：解出 $\eta$。** 把 (5) 与 (3) 代入 (34) 的系数：

$$
\frac{g\mu}{\pi\hbar^2k_za_\perp^2}=\frac{2a}{k_za_\perp^2},
\qquad
\frac{g\mu}{2\pi\hbar^2a_\perp}=\frac{a}{a_\perp},
$$

并注意 $\mathcal L=-\mathcal C+\bar{\mathcal L}$，把 (34) 写成

$$
\eta\left[1+\frac{2ia}{k_za_\perp^2}+\frac{a}{a_\perp}\big(\bar{\mathcal L}(\epsilon)-\mathcal C\big)\right]
=\frac{1}{\sqrt{\pi}a_\perp},
$$

故

$$
\eta=\frac{1}{\sqrt{\pi}a_\perp}\Big/
\left[1+\frac{2ia}{k_za_\perp^2}+\frac{a}{a_\perp}\big(\bar{\mathcal L}(\epsilon)-\mathcal C\big)\right].
\tag{39}
$$

**第二步：用 $\eta$ 表示 $f_{\rm even}$。** 由 (29)，开放道波函数在 $|z|\to\infty$ 的散射部分给出

$$
\psi_0(z)\xrightarrow{z\to+\infty}
e^{ik_zz}\Big[1-\frac{i\mu g}{\hbar^2k_z}\phi_0^*(0)\eta\Big].
$$

与 (17) 的渐近形式比较（横向部分 $\phi_0$ 提出后）：

$$
\boxed{\;
f_{\rm even}(k_z)=-\frac{i\mu g}{\hbar^2k_z}\phi_0^*(0)\,\eta
=-\frac{2i\sqrt{\pi}a}{k_za_\perp}\,\eta\;}
\tag{40}
$$

（第二步用了 $\phi_0(0)=1/(\sqrt{\pi}a_\perp)$ 与 $\mu g/\hbar^2=2\pi a$：$-i\frac{\mu g}{\hbar^2k_z}\frac{1}{\sqrt{\pi}a_\perp}=-\frac{2i\sqrt{\pi}a}{k_za_\perp}$。）

把 (39) 代入 (40)，分子分母同除以 $N=\frac{2ia}{k_za_\perp^2}$，得到**原文的核心结果（式 (11)）**：

$$
\boxed{\;
f_{\rm even}(k_z)=
-\frac{1}{\,1+ik_za_{\rm 1D}
-\dfrac{ik_za_\perp}{2}\,\bar{\mathcal L}\!\left(-\dfrac{k_z^2a_\perp^2}{4}\right)\,}
\;}
\tag{41}
$$

其中

$$
\boxed{\;
a_{\rm 1D}=-\frac{a_\perp^2}{2a}\Big(1-\mathcal C\frac{a}{a_\perp}\Big)
=-\frac{a_\perp^2}{2a}+\frac{\mathcal C}{2}a_\perp
\;}
\tag{42}
$$

（原文式 (12)）称为**一维散射长度**。

**代数细节（保证可复现）**：由 (39)(40) 有

$$
f_{\rm even}=-\frac{N}{1+N+\frac{a}{a_\perp}(\bar{\mathcal L}-\mathcal C)},
\qquad
N\equiv\frac{2ia}{k_za_\perp^2}.
$$

上下同除 $N$：

$$
f_{\rm even}=-\frac{1}{\frac1N+1+\frac{a}{a_\perp}(\bar{\mathcal L}-\mathcal C)\frac1N},
\qquad
\frac1N=\frac{k_za_\perp^2}{2ia}=-\frac{ik_za_\perp^2}{2a},
$$

$$
\frac{a}{a_\perp}\big(\bar{\mathcal L}-\mathcal C\big)\frac1N
=\frac{a}{a_\perp}\big(\bar{\mathcal L}-\mathcal C\big)\Big(-\frac{ik_za_\perp^2}{2a}\Big)
=-\frac{ik_za_\perp}{2}\big(\bar{\mathcal L}-\mathcal C\big).
$$

于是

$$
f_{\rm even}=-\frac{1}{1-\frac{ik_za_\perp^2}{2a}-\frac{ik_za_\perp}{2}\bar{\mathcal L}+\frac{ik_za_\perp\mathcal C}{2}}
=-\frac{1}{1+ik_z\Big(\underbrace{-\frac{a_\perp^2}{2a}+\frac{\mathcal C a_\perp}{2}}_{=a_{\rm 1D}}\Big)-\frac{ik_za_\perp}{2}\bar{\mathcal L}}\ ,
$$

即 (41)(42)。∎

**一维散射长度的另一种定义（原文）**：把偶道解写成 $\Psi\propto\sin[k_z|z|+\Delta(k_z)]\phi_{0,0}$，定义

$$
a_{\rm 1D}\equiv-\left.\frac{\partial\Delta}{\partial k_z}\right|_{k_z\to0^+}.
\tag{43}
$$

可以验证它与 (42) 一致：由 (41) 可得（见 §4.1 的推导思路）

$$
\tan\Delta=-k_za_{\rm 1D}+O(k_z^3)
\quad\Longrightarrow\quad
\Delta\simeq-k_za_{\rm 1D}
\quad\Longrightarrow\quad
-\partial_{k_z}\Delta\big|_{0}=a_{\rm 1D} .
$$

**重要物理量**：

- $\bar{\mathcal L}$ 项 $\sim O\big((k_za_\perp)^3\big)$（因为 $\bar{\mathcal L}=O(k_z^2a_\perp^2)$ 再乘 $k_za_\perp$），在低能下是小量；
- $a_{\rm 1D}$ 由 $a$、$a_\perp$ 完全确定：横向约束的全部低能效应都被"吸收"进一个参数；
- 当 $\mathcal C a/a_\perp=1$（即 $a=a_\perp/\mathcal C$）时 $a_{\rm 1D}=0$。由于 $g_{\rm 1D}\propto1/a_{\rm 1D}$（§4.2），这对应相互作用强度发散——**约束诱导共振（CIR）**。

---

## 4. 低能极限：有效一维 $\delta$ 势与约束诱导共振

### 4.1 一维 $\delta$ 势的严格解（推导）

先解一个独立的一维问题：质量为 $\mu$ 的粒子在一维 $\delta$ 势中散射。由上一节的结果知道，一维相互作用应该由 $a_{\rm 1D}$ 描述，因此我们设

$$
\boxed{\;
-\frac{\hbar^2}{2\mu}\psi''(z)+g_{\rm 1D}\delta(z)\psi(z)=E\psi(z),
\qquad E=\frac{\hbar^2k_z^2}{2\mu}\;}
\tag{44}
$$

（原文式 (13)），其中 $g_{\rm 1D}$ 待定。

**边界条件（推导）**：把方程在 $(-\eta,\eta)$ 上积分，取 $\eta\to0^+$：

$$
-\frac{\hbar^2}{2\mu}\big[\psi'(0^+)-\psi'(0^-)\big]+g_{\rm 1D}\psi(0)=0
\quad\Longrightarrow\quad
\boxed{\;\psi'(0^+)-\psi'(0^-)=\frac{2\mu g_{\rm 1D}}{\hbar^2}\,\psi(0)\;}
\tag{45}
$$

**求解**：取散射形式

$$
\psi(z)=e^{ik_zz}+f\,e^{ik_z|z|}.
$$

由 $z=0$ 处波函数连续：$\psi(0)=1+f$。求导：

$$
\psi'(0^+)=ik_z(1+f),\qquad
\psi'(0^-)=ik_z(1-f),
$$

代入 (45)：

$$
2ik_zf=\frac{2\mu g_{\rm 1D}}{\hbar^2}(1+f)
\ \Longrightarrow\
f=\frac{\mu g_{\rm 1D}/\hbar^2}{ik_z-\mu g_{\rm 1D}/\hbar^2}.
$$

令 $\alpha\equiv\mu g_{\rm 1D}/\hbar^2$，则

$$
f=\frac{\alpha}{ik_z-\alpha}=-\frac{1}{1-ik_z/\alpha}
=-\frac{1}{1+ik_z(-\frac1\alpha)} .
$$

与低能一维振幅的标准形式 $f=-1/(1+ik_za_{\rm 1D})$ 比较，得到

$$
\boxed{\;
f^{\delta}(k_z)=-\frac{1}{1+ik_za_{\rm 1D}},
\qquad
a_{\rm 1D}=-\frac{\hbar^2}{\mu g_{\rm 1D}}
\ \Longleftrightarrow\
g_{\rm 1D}=-\frac{\hbar^2}{\mu a_{\rm 1D}} .
\;}
\tag{46}
$$

**一致性与相位的检验**：把偶道解写成 $\psi\propto\sin(k_z|z|+\Delta)$。对 (46) 的解可验证
$\tan\Delta=-k_za_{\rm 1D}$，故 $\Delta\simeq-k_za_{\rm 1D}$，与 (43) 的定义 $a_{\rm 1D}=-\partial_{k_z}\Delta|_{0}$ 一致。∎

**透射系数（推导）**：由渐近形式 (17)，透射波的振幅是 $1+f_{\rm even}+f_{\rm odd}=1+f$（在一维中 $z\to+\infty$ 的波）。因此

$$
\boxed{\;
\mathcal T=\big|1+f_{\rm even}+f_{\rm odd}\big|^2
=\big|1+f\big|^2
=\frac{(k_za_{\rm 1D})^2}{1+(k_za_{\rm 1D})^2}
\ \xrightarrow{k_z\to0}\ 0 .\;}
\tag{47}
$$

**结论**：一维 $\delta$ 势在零能极限下**全反射**（只要 $a_{\rm 1D}\neq0$）；这是后文 Tonks 气体（硬核）的根源。

### 4.2 有效一维相互作用强度 $g_{\rm 1D}$（推导）

把精确振幅 (41) 在低能展开。因为 $\bar{\mathcal L}=O((k_za_\perp)^2)$，所以

$$
f_{\rm even}(k_z)=-\frac{1}{1+ik_za_{\rm 1D}+O\big((k_za_\perp)^3\big)}
\ \xrightarrow{k_za_\perp\ll1}\
-\frac{1}{1+ik_za_{\rm 1D}}=f^\delta(k_z).
$$

即**低能下精确理论等价于一个一维 $\delta$ 势** (44)。对照 (46) 与 (42)：

$$
\boxed{\;
g_{\rm 1D}=-\frac{\hbar^2}{\mu a_{\rm 1D}}
=\frac{2\hbar^2a}{\mu a_\perp^2}\,
\frac{1}{1-\mathcal C a/a_\perp}
=g\,|\phi_{0,0}(0)|^2\,
\frac{1}{1-\mathcal C a/a_\perp}\;}
\tag{48}
$$

（原文式 (14)）。**推导最后一步**：

$$
-\frac{\hbar^2}{\mu a_{\rm 1D}}
=-\frac{\hbar^2}{\mu}\cdot\frac{1}{-\frac{a_\perp^2}{2a}\big(1-\mathcal C a/a_\perp\big)}
=\frac{2\hbar^2a}{\mu a_\perp^2}\frac{1}{1-\mathcal C a/a_\perp},
$$

并用 (5)(12)：$g|\phi_{0,0}(0)|^2=\frac{2\pi\hbar^2a}{\mu}\cdot\frac{1}{\pi a_\perp^2}=\frac{2\hbar^2a}{\mu a_\perp^2}$。∎

**物理意义（难点解说）**：

- $g|\phi_{0,0}(0)|^2$ 是"把三维耦合按横向基态在原点的概率权重摊开"的朴素结果；
- 分母 $(1-\mathcal C a/a_\perp)^{-1}$ 是**横向虚激发**（式 (31) 第 (III) 项）对相互作用的**重正化**：常数 $\mathcal C$ 是无穷多闭道求和（$\zeta$ 函数）产生的；
- $a\to0$（弱相互作用）时 $g_{\rm 1D}\approx 2\hbar^2a/(\mu a_\perp^2)\to0$，符合直觉；
- $a\to a_\perp/\mathcal C$ 时 $g_{\rm 1D}\to\infty$：约束诱导共振。

### 4.3 约束诱导共振（CIR）

由 (48) 分母为零的条件

$$
\boxed{\;1-\mathcal C\frac{a}{a_\perp}=0
\ \Longleftrightarrow\
a=\frac{a_\perp}{\mathcal C}\approx0.685\,a_\perp\;}
\tag{49}
$$

给出**约束诱导共振（confinement-induced resonance, CIR）**：即使真实的三维散射长度 $a$ 有限（没有三维共振），横向约束本身也会在 $a=a_\perp/\mathcal C$ 处产生一维共振（$g_{\rm 1D}\to\pm\infty$、$a_{\rm 1D}=0$）。

**物理解释**：碰撞过程中被虚激发的横向模式（闭道云）相当于把相互作用"改造"了。当横向束缚能与三维散射的准束缚态能量匹配时，出现共振。CIR 的位置由约束几何（$a_\perp$）决定，因此实验上可以通过改变陷阱频率或利用 Feshbach 共振调节 $a$ 来扫描 CIR。

**符号与图 1**：

- $a<a_\perp/\mathcal C$（且 $a>0$）：$g_{\rm 1D}>0$（排斥），$a_{\rm 1D}<0$；
- $a>a_\perp/\mathcal C$ 或 $a<0$：$g_{\rm 1D}<0$（吸引）。
- 原文图 1 比较了精确透射系数 (47)（实线）与 $\delta$ 势近似（虚线）随入射动量 $k_z$ 的变化，参数取 $a_\perp/a=\pm10$ 等。当 $|a|$ 小（弱相互作用）时两者几乎重合；在低能端（$k_za_\perp\ll1$）也都趋于全反射。**图 1 的作用**：验证低能展开 (48) 的适用范围——即"何时可以用有效一维 $\delta$ 势"。

### 4.4 与硬球散射的类比

原文指出：有效势 (46)(48) 可以证明能**重现半径为 $a$ 的硬球在横向陷阱中**的低能散射性质。理解这一点只需回忆：

- 在自由空间中，Huang 赝势 (4) 正是用来"零程地"重现半径为 $a$ 的硬球的低能散射（§2.2）；
- 在横向陷阱中，类比成立：低能散射由 $a_{\rm 1D}$（等价地 $g_{\rm 1D}$）完全描述；
- 因此与自由空间完全类似（原文引 Huang 教材），(48) 也可以推广到**多体问题**：只要所有原子处于横向基态、气体足够稀薄（三体碰撞可忽略），准一维多体系统的相互作用即可用 $U_{\rm 1D}(z)=g_{\rm 1D}\delta(z)$ 描述。这直接导向下一节的 Tonks 气体。

---

## 5. 多体应用：Tonks 气体与玻色–费米对偶

### 5.1 强耦合极限：全反射与不可穿透

由 (47)，当

$$
\boxed{\;k_z|a_{\rm 1D}|\ll1\;}
\tag{50}
$$

（原文式 (15)）时透射系数趋于零，即**全反射**。在低能散射极限下，系统的两体散射等价于**不可穿透的硬核**（全反射）相互作用。

**区别（易混点）**：

- CIR 条件 (49) 说的是"横向约束使 $g_{\rm 1D}$ 发散"，即 $a_{\rm 1D}=0$；
- 全反射条件 (50) 说的是"能量足够低"，即使 $g_{\rm 1D}$ 有限，$k_z\to0$ 时透射也趋于零。
- 对**多体**气体，$k_z$ 的最大值由费米动量给出（§5.5），于是 (50) 变成对密度/原子数的限制。

### 5.2 Girardeau 玻色–费米映射（推导）

**多体模型**：N 个质量为 $m$ 的玻色子（无自旋）在一维盒中，通过 (48) 的有效相互作用相互作用：

$$
H=-\frac{\hbar^2}{2m}\sum_{i=1}^{N}\frac{\partial^2}{\partial z_i^2}
+g_{\rm 1D}\sum_{i<j}\delta(z_i-z_j),
\qquad g_{\rm 1D}\to+\infty\ (\text{硬核极限}),
\tag{51}
$$

波函数对任意坐标交换**对称**。硬核极限下，两体接触边界条件 (45) 变为

$$
\boxed{\;\Psi(z_1,\ldots,z_N)\big|_{z_i=z_j}=0\;}
\tag{52}
$$

**映射构造**：取 N 个无自旋无相互作用费米子的基态 $\Psi^f$（反对称），定义

$$
\boxed{\;\Psi^{b}=\big|\Psi^{f}\big|\;}
\tag{53}
$$

（原文式 (16)）。逐条验证：

1. **对称性**：$|\Psi^f|$ 在交换下不变 ✓。
2. **硬核条件**：$\Psi^f$ 反对称，任意两坐标重合时 $\Psi^f=0$，故 $|\Psi^f|=0$ ✓ 满足 (52)。
3. **运动方程**：把构型空间按坐标大小分成 $N!$ 个区域（例如 $z_1<\cdots<z_N$ 及其置换）。在每个区域内部 $|\Psi^f|=\pm\Psi^f$，满足自由方程；$\delta$ 项只在区域边界（重合点）出现。故 $\Psi^b$ 在每个区域内部满足自由薛定谔方程，其能量与 $\Psi^f$ 相同。
4. **边界连接**：无限大耦合的接触条件就是 $\Psi=0$（由 (52)），而 $|\Psi^f|$ 恰好在边界为零。导数的跳变由接触条件控制（对无穷强耦合，跳变条件无矛盾）。
5. **基态**：可以证明每个玻色本征态都对应一个费米本征态（取绝对值），映射保持能量，因此基态能量等于理想费米气体的基态能量。

**为什么只有一维**（原文引文 [7]）：在 $d\ge2$ 中，两粒子"重合"的超曲面不能把构型空间分割成互不相连的排序区域，粒子可以"绕过"彼此；$|\Psi^f|$ 无法拼成一个在全局满足自由方程的对称波函数。一维的特殊拓扑是 Girardeau 映射成立的关键。

**物理后果（玻色–费米对偶/fermionization）**：Tonks 气体的**能量与热力学**和理想（无自旋）费米气体相同；但**关联函数**与**动量分布**不同（因为它依赖于取绝对值这一操作）。

### 5.3 基态与费米海（推导）

取 N 个费米子，长度 $L$ 的盒子，周期边界条件。单粒子平面波态

$$
\langle z|j\rangle=\frac{1}{\sqrt L}e^{ik_jz},
\qquad
k_j=\frac{2\pi}{L}j,\quad j\in\mathbb Z .
$$

基态 Slater 行列式（原文式 (17)）：

$$
\boxed{\;
\Psi^{f}(z_1,\ldots,z_N)
=\frac{1}{\sqrt{N!\,L^{N}}}\det\!\big(e^{ik_jz_{j'}}\big)_{j,j'=1}^{N}
\;}
\tag{54}
$$

**归一化验证**：由轨道正交性 $\int_0^Ldz\,e^{i(k_j-k_{j'})z}=L\delta_{j,j'}$，标准行列式积分给出 $\int d^Nz|\Psi^f|^2=1$。∎

**费米海**：把最低的 N 个态填满，取对称分布 $j=-(N-1)/2,\ldots,(N-1)/2$（$N$ 为奇数；偶数时差 $O(1)$）。于是

$$
\boxed{\;k_{\rm Fermi}=\frac{\pi(N-1)}{L}\;}
\tag{55}
$$

（原文式 (18)）。基态能量

$$
E^b=E^f=\sum_{j\ \rm occ}\frac{\hbar^2k_j^2}{2m}
$$

就是 Tonks 气体的基态能量。

### 5.4 动量分布（引用已知结果与解释）

**定义**：一体动量分布 $w(k_z)$ 由单体密度矩阵 $\rho_1(z-z')=\langle\Psi^\dagger(z)\Psi(z')\rangle$ 的 Fourier 变换给出；原文归一化为 $\int(dk_z/k_{\rm Fermi})w(k_z)=1$。

**原文图 2 所展示的已知结果**（来自引文 [14]：Vaidya–Tracy 1979，修正版 Jimbo 等 1980。这一分布不是本文推导的，而是在文中引用并作图）：

1. **小动量端有 $1/\sqrt{k_z}$ 发散**：

   $$
   \boxed{\;w(k_z)\simeq\rho_\infty\sqrt{\frac{k_{\rm Fermi}}{2\pi k_z}}\;}
   \tag{56}
   $$

   $$
   \rho_\infty=\pi e^{1/2}2^{-1/3}A^{-6}=0.92418\ldots,
   \qquad A=1.2824\ldots\ \text{(Glaisher 常数)}.
   $$
2. **大动量端**为幂律尾（$\propto1/k_z^4$，来自波函数在接触点的尖点/节点）。
3. 与理想费米气体（$T=0$ 时为阶跃函数）对比：Tonks 气体的分布被"抹平"，并且在原点出现 $1/\sqrt{k_z}$ 峰。

**物理解释（难点解说）**：

- 费米阶跃的消失意味着"费米面"在动量分布中不再清晰；这是强关联效应；
- $1/\sqrt{k_z}$ 峰是 Tonks 气体（准长程关联/fermionization）的**实验判据**：只要在动量分布中看到原点处的 $1/\sqrt{k_z}$ 峰，就说明系统进入了不可穿透区；
- 由于这一发散发生在**热力学极限**，对有限 $(N,L)$ 而言，在 $k_z\lesssim2\pi/L\sim2k_{\rm Fermi}/N$ 的区域分布曲线失效（热力学极限在该小动量区不成立）；
- 大动量尾的物理来源是接触点的节点（波函数在 $z_i=z_j$ 处线性消失，导致密度矩阵在短距离处有 $|z|$ 型行为，其 Fourier 变换为 $1/k_z^4$）。

### 5.5 原子数上限与实验参数（推导）

**推导 $N^*$**：多体气体中，任意一对原子的**相对动量**为

$$
p_{\rm rel}=\mu(v_2-v_1)=\frac m2(v_2-v_1)
\ \Longrightarrow\
k_{\max}=\frac{\mu}{m}\max_{j,j'}|k_j-k_{j'}|
=\frac12\cdot2k_{\rm Fermi}=k_{\rm Fermi}.
$$

（原文式 (19) 前的论证。）把全反射条件 (50) 用到**所有**原子对，即 $k_{\max}|a_{\rm 1D}|\ll1$，并用 (55)：

$$
\frac{\pi(N-1)}{L}|a_{\rm 1D}|\ll1
\ \Longrightarrow\
\boxed{\;N\ll N^*=\frac{L}{\pi|a_{\rm 1D}|}\;}
\tag{57}
$$

（原文式 (19)）。由于（当 $a\ll a_\perp$ 时）$|a_{\rm 1D}|\approx a_\perp^2/(2a)\propto1/\omega_\perp$，所以

$$
\boxed{\;N^*\propto L\,\omega_\perp\;}
$$

**约束越强、盒子越长，可容纳的 Tonks 气体原子越多**。

**数值检验（代入 (3)(42)(57)）**：取 $\omega_\perp=2\pi\times10^4\,$Hz、$L=3\,$cm：

| 原子        | $a$                  | $a_\perp$                                   | $\|a_{\rm 1D}\|$   | $N^*$         |
| ----------- | ---------------------- | --------------------------------------------- | -------------------- | --------------- |
| $^{87}$Rb | $+110\,a_{\rm Bohr}$ | $2.9\times10^3a_{\rm Bohr}\approx1.5\,\mu$m | $\approx1.9\,\mu$m | $5\times10^3$ |
| Na          | $+52\,a_{\rm Bohr}$  | $5.6\times10^3a_{\rm Bohr}\approx3.0\,\mu$m | $\approx16\,\mu$m  | $6\times10^2$ |

**计算示例（Rb）**：$a_\perp=\sqrt{2\hbar/(m\omega_\perp)}=2.9\times10^3a_{\rm Bohr}$；
$1-\mathcal C a/a_\perp=1-1.460\times110/2900=0.945$；
$|a_{\rm 1D}|=\frac{a_\perp^2}{2a}(1-\mathcal C a/a_\perp)\approx\frac{(2900)^2}{220}\times0.945\,a_{\rm Bohr}\approx3.6\times10^4a_{\rm Bohr}\approx1.9\,\mu$m；
$N^*=L/(\pi|a_{\rm 1D}|)=0.03\,\text{m}/(\pi\times1.9\,\mu\text{m})\approx5\times10^3$ ✓（与原文一致；Na 同理得 $6\times10^2$）。

### 5.6 实验实现条件（原文结论）

综合全文，实现 Tonks 气体的实验条件为：

1. **高度拉长的陷阱（波导）**：横向强束缚 $\omega_\perp$ 大，纵向自由；
2. **低温**：$k_BT\ll\hbar\omega_\perp$，保证横向冻结在基态；
3. **低线性密度**：$\rho a\ll(2\pi)^{-1}(a/a_\perp)^2$（$a>0$）——这是 $k_F|a_{\rm 1D}|\ll1$ 的数量级条件，保证费米海内所有碰撞都在全反射区；
4. **正的三维散射长度** $a>0$ 且 $a<a_\perp/\mathcal C$（对应排斥的有效一维相互作用 $g_{\rm 1D}>0$）。

实验方案：**原子波导 + 纵向势垒**（长度 $L$ 的盒）；判据：动量分布中的 $1/\sqrt{k_z}$ 峰（图 2）。低能条件 (50) 同时给出原子数上限 (57)。

---

## 6. 适用范围、总结与后记

### 6.1 赝势近似的适用条件（原文讨论）

**(i) 振幅 (41) 的适用范围**：推导 (41) 时只用了"相互作用是零程赝势"与"横向模式仅 $m=0$、只有基态开放"，没有对横向约束的强弱或 $a$ 的大小做额外假设。因此 (41) 对**任意强度的横向约束**成立。

**(ii) 紧束缚极限 $a_\perp\ll|a|$ 的前提**：这个极限只有在真实三维势处于**零能共振**附近时才有意义，即

$$
|a|\gg r_0\qquad(r_0:\ \text{真实势的有效力程}).
\tag{58}
$$

**推导（条件链）**：

1. 赝势近似 (4) 只在 $k\ll r_0^{-1}$ 时有效（超出这个动量，真实势的形状开始起作用；见 §2.2 的适用条件）。
2. 在横向陷阱中，碰撞过程中涉及的横向动量尺度为 $k_\perp\sim1/a_\perp$，因此要求

   $$
   \frac{1}{a_\perp}\ll\frac1{r_0}\quad\Longleftrightarrow\quad a_\perp\gg r_0 .
   $$
3. 同时，紧束缚（横向冻结）要求 $a_\perp\ll|a|$。
4. 两者相容要求 $r_0\ll a_\perp\ll|a|$，即

   $$
   \boxed{\;|a|\gg r_0\;}
   $$

   这正是零能共振条件。满足它的典型原子：Cs、$^{85}$Rb、$^{39}$K（原文引文 [11]）。
5. 若 $|a|$ 不够大，则不能同时满足"赝势有效"与"横向冻结"，本文的零程模型失效。

**(iii) 其他隐含条件**：

- 温度足够低（只有横向基态）：$k_BT\ll\hbar\omega_\perp$；
- 稀薄气体（忽略三体碰撞）：这是把两体有效势用于多体问题的前提（§4.4）；
- 无自旋/自旋极化（或自旋自由度冻结），本文的赝势是自旋无关的。

### 6.2 结果总结

| 结果                | 公式                                                                                                       | 本材料编号 | 原文编号 |
| ------------------- | ---------------------------------------------------------------------------------------------------------- | ---------- | -------- |
| 精确一维散射振幅    | $f_{\rm even}=-\big[1+ik_za_{\rm 1D}-\frac{ik_za_\perp}{2}\bar{\mathcal L}(-k_z^2a_\perp^2/4)\big]^{-1}$ | (41)       | (11)     |
| 一维散射长度        | $a_{\rm 1D}=-\frac{a_\perp^2}{2a}(1-\mathcal C a/a_\perp)$                                               | (42)       | (12)     |
| 有效一维势          | $U_{\rm 1D}=g_{\rm 1D}\delta(z)$                                                                         | (44)       | (13)     |
| 有效耦合            | $g_{\rm 1D}=-\frac{\hbar^2}{\mu a_{\rm 1D}}=g|\phi_{0,0}(0)|^2(1-\mathcal C a/a_\perp)^{-1}$             | (48)       | (14)     |
| CIR 位置            | $a=a_\perp/\mathcal C$，$\mathcal C=1.4603\ldots=-\zeta(1/2)$                                          | (49)       | —       |
| 全反射/不可穿透条件 | $k_z|a_{\rm 1D}|\ll1$                                                                                    | (50)       | (15)     |
| 玻色–费米映射      | $\Psi^b=|\Psi^f|$                                                                                        | (53)       | (16)     |
| 费米海              | $k_{\rm Fermi}=\pi(N-1)/L$                                                                               | (55)       | (18)     |
| 动量分布小$k$ 峰  | $w\simeq\rho_\infty\sqrt{k_{\rm Fermi}/(2\pi k_z)}$                                                      | (56)       | 图 2     |
| 原子数上限          | $N\ll N^*=L/(\pi|a_{\rm 1D}|)\propto L\omega_\perp$                                                      | (57)       | (19)     |

**一句话总结**：横向简谐约束把三维低能散射约化为一个一维 $\delta$ 势问题；其耦合常数由 $a$ 和 $a_\perp$ 通过 $a_{\rm 1D}$ 决定，且可以被约束诱导共振放大到无穷；在低能（低密度）极限下，准一维玻色气体成为 Tonks 气体——一个与理想费米气体一一对应的强关联一维系统。

### 6.3 后记：原文之后的实验进展（超出原文范围，仅供参考）

本文 1998 年发表后，其预言在冷原子实验中得到了检验：

- **约束诱导共振**：在准一维 Cs 气体中被观测到（Moritz 等，Phys. Rev. Lett. **91**, 250402 (2003)）。
- **Tonks 气体**：在光学晶格/波导中的准一维玻色气体中实现（Paredes 等，Nature **429**, 277 (2004)；Kinoshita 等，Science **305**, 1125 (2004)）。
- 动量分布的 $1/\sqrt{k_z}$ 特征、玻色–费米对偶的热力学后果等随后也被实验研究。

（这些内容不属于原论文的推导范围，列出仅供建立物理背景。）

---

## 附录 A：二维谐振子（本征函数与 $|\phi_{n,0}(0)|^2$）

### A.1 笛卡尔坐标解

由 (9)：

$$
\hat H_\perp=-\frac{\hbar^2}{2\mu}\Big(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\Big)
+\frac{\mu\omega_\perp^2}{2}(x^2+y^2)
=\hat h(x)+\hat h(y),
$$

$$
\hat h(s)=-\frac{\hbar^2}{2\mu}\frac{d^2}{ds^2}+\frac{\mu\omega_\perp^2}{2}s^2 .
$$

一维谐振子本征值 $\hbar\omega_\perp(n_s+\tfrac12)$，因此

$$
E_{n_x,n_y}=\hbar\omega_\perp(n_x+n_y+1)\equiv\hbar\omega_\perp(n+1),
\qquad n=0,1,2,\ldots
\tag{A.1}
$$

与 (10) 一致。基态（$n_x=n_y=0$）：

$$
\phi_{0,0}(x,y)=\frac{1}{\sqrt{\pi}a_\perp}e^{-(x^2+y^2)/(2a_\perp^2)},
\qquad a_\perp=\sqrt{\frac{\hbar}{\mu\omega_\perp}},
\tag{A.2}
$$

即 (11)。

### A.2 极坐标、径向方程与 $m=0$ 解（推导）

用极坐标 $(\rho,\varphi)$：

$$
\hat H_\perp=-\frac{\hbar^2}{2\mu}\Big[\frac1\rho\frac{\partial}{\partial \rho}\Big(\rho\frac{\partial}{\partial\rho}\Big)
+\frac{1}{\rho^2}\frac{\partial^2}{\partial\varphi^2}\Big]
+\frac{\mu\omega_\perp^2}{2}\rho^2 .
$$

设 $\phi=R(\rho)e^{im\varphi}/\sqrt{2\pi}$，径向方程为

$$
-\frac{\hbar^2}{2\mu}\Big[R''+\frac1\rho R'-\frac{m^2}{\rho^2}R\Big]
+\frac{\mu\omega_\perp^2}{2}\rho^2R=ER .
\tag{A.3}
$$

代入 $E=\hbar\omega_\perp(n+1)$，并令 $x=\rho/a_\perp$：

$$
R_{xx}+\frac1xR_x-\frac{m^2}{x^2}R
+\big[2(n+1)-x^2\big]R=0 .
\tag{A.4}
$$

**$m=0$ 情形（本文用到的唯一情形）**：方程变为

$$
R_{xx}+\frac1xR_x+\big[2(n+1)-x^2\big]R=0 .
\tag{A.5}
$$

作 $R=e^{-x^2/2}f(x^2)$ 可化为拉盖尔方程（见下），解为

$$
\boxed{\;
R_{n,0}(\rho)=\sqrt{\frac{2}{a_\perp^2}}\,
e^{-\rho^2/(2a_\perp^2)}\,
L_{n/2}\!\Big(\frac{\rho^2}{a_\perp^2}\Big),
\qquad n=0,2,4,\ldots\;}
\tag{A.6}
$$

其中 $L_k$ 是拉盖尔多项式，满足

$$
t\,L_k''+(1-t)L_k'+kL_k=0,\qquad L_k(0)=1,\qquad
\int_0^\infty L_k(t)L_j(t)e^{-t}dt=\delta_{kj}.
\tag{A.7}
$$

**验证（$k=0$ 与 $k=1$）**：

- $k=0$：$R=e^{-x^2/2}$。直接代入 (A.5)（$n=0$）：$R_{xx}+R_x/x+(2-x^2)R=(x^2-1)e^{-x^2/2}-e^{-x^2/2}+(2-x^2)e^{-x^2/2}=0$ ✓。
- $k=1$（$n=2$）：$R=e^{-x^2/2}(1-x^2)$。代入 (A.5)（$n=2$，即 $[6-x^2]$）：由 $g=1-x^2$、$g'=-2x$、$g''=-2$，得 $R''=e^{-x^2/2}(g''-2xg'+(x^2-1)g)=e^{-x^2/2}(-3+6x^2-x^4)$；$R'/x=e^{-x^2/2}(-3+x^2)$。相加：$e^{-x^2/2}(-6+7x^2-x^4)$；再加上 $(6-x^2)R=e^{-x^2/2}(6-7x^2+x^4)$，总和为零 ✓。

**归一化**：令 $t=\rho^2/a_\perp^2$（$\rho d\rho=\frac{a_\perp^2}{2}dt$）：

$$
\int_0^\infty |R_{n,0}|^2\rho\,d\rho
=\frac{2}{a_\perp^2}\cdot\frac{a_\perp^2}{2}
\int_0^\infty L_k(t)^2e^{-t}dt=1
\qquad(k=n/2),
$$

用了 (A.7) 的正交归一。于是全波函数

$$
\phi_{n,0}(\rho)=\frac{1}{\sqrt{2\pi}}R_{n,0}(\rho),
\qquad
\int|\phi_{n,0}|^2d^2\rho=1 .
\tag{A.8}
$$

### A.3 $|\phi_{n,0}(0)|^2$

由 (A.6)(A.8) 与 $L_k(0)=1$：

$$
\phi_{n,0}(0)=\frac{1}{\sqrt{2\pi}}\sqrt{\frac{2}{a_\perp^2}}=\frac{1}{\sqrt{\pi}a_\perp},
\qquad
\boxed{\;|\phi_{n,0}(0)|^2=\frac{1}{\pi a_\perp^2}\quad(\forall n\ \text{偶})\;}
\tag{A.9}
$$

即 (13)。

---

## 附录 B：一维 Green 函数的详细推导

### B.1 出射波 Green 函数 $G_0$

要求 $(\partial_z^2+k^2)G_0(z)=\delta(z)$ 且 $|z|\to\infty$ 为出射波。Fourier 表示：

$$
G_0(z)=\int_{-\infty}^{\infty}\frac{dq}{2\pi}\frac{e^{iqz}}{k^2-q^2+i0^+}.
$$

被积函数在 $q=\pm(k+i0^+)$ 有极点。对 $z>0$ 把围道在上半平面闭合（$e^{iqz}$ 在上半平面衰减），围住的极点是 $q=+k+i0^+$。把分母分解：

$$
k^2-q^2=(k-q)(k+q),\qquad
q\simeq k\ \text{处}\ k+q\simeq2k,
$$

故留数为

$$
\underset{q=k}{\mathrm{Res}}\,\frac{e^{iqz}}{k^2-q^2}
=\frac{e^{ikz}}{2k}\cdot(-1)=-\frac{e^{ikz}}{2k},
$$

于是 $G_0(z>0)=2\pi i\cdot\frac{1}{2\pi}\left(-\frac{e^{ikz}}{2k}\right)=-\frac{i}{2k}e^{ikz}$；对 $z<0$ 同理（换成 $e^{-ikz}$）。合起来

$$
G_0(z)=-\frac{i}{2k}e^{ik|z|}.
\tag{B.1}
$$

**验证**：$G_0'(z)=-\frac{i}{2k}\cdot ik\,\mathrm{sign}(z)e^{ik|z|}=\frac12\mathrm{sign}(z)e^{ik|z|}$，故 $G_0'(0^+)-G_0'(0^-)=1$；且对 $z\neq0$，$(\partial_z^2+k^2)G_0=0$。于是 $(\partial_z^2+k^2)G_0=\delta(z)$ ✓。

### B.2 衰减 Green 函数 $G_s$

要求 $(\partial_z^2-\kappa^2)G_s=\delta(z)$、$|z|\to\infty$ 衰减：

$$
G_s(z)=-\frac{1}{2\kappa}e^{-\kappa|z|} .
\tag{B.2}
$$

**验证**：$G_s'=\frac12\mathrm{sign}(z)e^{-\kappa|z|}$，跳变 $=1$；对 $z\neq0$，$(\partial_z^2-\kappa^2)G_s=0$ ✓。

### B.3 用法

对

$$
(\partial_z^2-\lambda^2)\psi(z)=S\delta(z),
\qquad
\lambda=\begin{cases}ik_z& (\text{开放道})\\ \kappa_s&(\text{闭合道})\end{cases},
$$

的解为 $\psi(z)=S\,G(z)$，其中 $G$ 分别取 (B.1)(B.2)（这正是 §3.3 中 (29)(30) 的来源）。

---

## 附录 C：$\Lambda$ 函数、$\zeta$ 函数与常数 $\mathcal C$（Mellin 变换推导）

### C.1 定义与 Mellin 变换

$$
\Lambda[\xi,\epsilon]=\sum_{s=1}^{\infty}\frac{e^{-\sqrt{s+\epsilon}\,\xi}}{\sqrt{s+\epsilon}},
\qquad \xi>0 .
\tag{C.1}
$$

对 $\xi$ 作 Mellin 变换：

$$
M(\nu)\equiv\int_0^\infty d\xi\,\xi^{\nu-1}\Lambda[\xi,\epsilon]
=\sum_{s=1}^{\infty}\frac{1}{\sqrt{s+\epsilon}}\int_0^\infty d\xi\,\xi^{\nu-1}e^{-\sqrt{s+\epsilon}\xi}
=\Gamma(\nu)\sum_{s=1}^{\infty}(s+\epsilon)^{-(\nu+1)/2},
$$

$$
\boxed{\;M(\nu)=\Gamma(\nu)\,\zeta\!\Big(\frac{\nu+1}{2},\,1+\epsilon\Big)\;}
\tag{C.2}
$$

其中 $\zeta(\alpha,q)=\sum_{n=0}^{\infty}(n+q)^{-\alpha}$ 是 **Hurwitz $\zeta$ 函数**（这里 $s-1\mapsto n$，$s+\epsilon=n+1+\epsilon$）。逆变换：

$$
\Lambda[\xi,\epsilon]=\frac{1}{2\pi i}\int_{c-i\infty}^{c+i\infty}d\nu\,M(\nu)\,\xi^{-\nu},
\qquad c>1 .
\tag{C.3}
$$

### C.2 小 $\xi$ 展开与 $\mathcal L(\epsilon)$

把围道向左闭合（$\xi<1$ 时 $\xi^{-\nu}$ 在左半平面衰减）。极点为：

1. **$\nu=1$**：$\zeta\big(\frac{\nu+1}{2},1+\epsilon\big)$ 的单极点。记 $w=(\nu+1)/2$，则 $w\to1$ 时 $\zeta(w,q)\simeq\frac{1}{w-1}=\frac{2}{\nu-1}$，故留数 $=\Gamma(1)\cdot2\cdot\xi^{-1}=\frac{2}{\xi}$。
2. **$\nu=0$**：$\Gamma(\nu)\simeq\frac1\nu$ 的单极点，留数 $=\zeta\big(\frac12,1+\epsilon\big)\xi^{0}=\zeta\big(\frac12,1+\epsilon\big)$。
3. **$\nu=-1,-2,\ldots$**：给出 $\mathcal O(\xi),\mathcal O(\xi^2),\ldots$ 的正规项。

于是

$$
\boxed{\;
\Lambda[\xi,\epsilon]=\frac{2}{\xi}+\mathcal L(\epsilon)+\mathcal O(\xi),
\qquad
\mathcal L(\epsilon)=\zeta\!\Big(\frac12,\,1+\epsilon\Big)\;}
\tag{C.4}
$$

**与原文对照**：原文写 $\mathcal L(\epsilon)=-\mathcal C+\bar{\mathcal L}(\epsilon)$，其中 $\mathcal C=-\zeta(\frac12)$。下面证明两者一致。

### C.3 $\epsilon$ 展开：$\mathcal C$ 与 $\bar{\mathcal L}$

Hurwitz $\zeta$ 对第二变量的导数公式：

$$
\frac{\partial^k}{\partial\epsilon^k}\zeta(\alpha,1+\epsilon)
=(-1)^k(\alpha)_k\,\zeta(\alpha+k,1+\epsilon),
\qquad
(\alpha)_k=\alpha(\alpha+1)\cdots(\alpha+k-1).
\tag{C.5}
$$

在 $\epsilon=0$ 处展开 $\mathcal L(\epsilon)=\zeta(\frac12,1+\epsilon)$：

$$
\mathcal L(\epsilon)=\sum_{k=0}^{\infty}\frac{\epsilon^k}{k!}\,(-1)^k\Big(\frac12\Big)_k\zeta\Big(\frac12+k\Big),
\qquad
\Big(\frac12\Big)_k=\frac{(2k-1)!!}{2^k} .
\tag{C.6}
$$

- $k=0$ 项：$\zeta(\frac12)=-1.4603545\ldots$。
- $k\ge1$ 项：与 (38) 完全一致。

因此

$$
\boxed{\;
\mathcal C=-\zeta\!\Big(\frac12\Big)=1.4603\ldots,
\qquad
\bar{\mathcal L}(\epsilon)=\sum_{k=1}^{\infty}(-1)^k
\frac{\zeta\big(\frac{1+2k}{2}\big)(2k-1)!!}{2^kk!}\epsilon^k\;}
\tag{C.7}
$$

**与原文定义 (37) 的等价性**：对 $0<s<1$，Riemann $\zeta$ 函数的解析延拓给出

$$
\zeta\!\Big(\frac12\Big)
=\lim_{S\to\infty}\left[\sum_{s'=1}^{S}\frac{1}{\sqrt{s'}}-\int_0^{S}\frac{ds'}{\sqrt{s'}}\right]
=\lim_{S\to\infty}\left[\sum_{s'=1}^{S}\frac{1}{\sqrt{s'}}-2\sqrt S\right],
$$

故

$$
\mathcal C=\lim_{S\to\infty}\left[2\sqrt S-\sum_{s'=1}^S\frac1{\sqrt{s'}}\right]=-\zeta\!\Big(\frac12\Big)
\tag{C.8}
$$

与原文式 (9) 相同（原文用 $\int_0^s$ 与 $\sum_{s'=1}^s$ 的差）。∎

### C.4 常用数值

$$
\zeta(1/2)=-1.4603545088\ldots,\quad
\zeta(3/2)=2.612375\ldots,\quad
\zeta(5/2)=1.341487\ldots .
$$

---

## 附录 D：符号表

| 符号                                                                                              | 含义                     | 量纲                         |
| ------------------------------------------------------------------------------------------------- | ------------------------ | ---------------------------- |
| $m$                                                                                             | 原子质量                 | 质量                         |
| $\mu=m/2$                                                                                       | 两体约化质量             | 质量                         |
| $\omega_\perp$                                                                                  | 横向谐振子频率           | 1/时间                       |
| $a_\perp=\sqrt{\hbar/(\mu\omega_\perp)}$                                                        | 横向基态尺度（振子长度） | 长度                         |
| $\rho=\sqrt{x^2+y^2}$                                                                           | 横向坐标                 | 长度                         |
| $z$                                                                                             | 纵向坐标                 | 长度                         |
| $k_z$ | 纵向相对动量（$\hbar k_z$ 为动量）                                                    | 1/长度                   |                              |
| $E=\hbar^2k_z^2/(2\mu)+\hbar\omega_\perp$                                                       | 相对运动总能量           | 能量                         |
| $\phi_{n,m}$                                                                                    | 二维谐振子本征态         | $1/$长度                   |
| $V(\mathbf r)$、$r_0$                                                                         | 真实两体势与其有效力程   | 能量、长度                   |
| $a$                                                | 三维$s$ 波散射长度                       | 长度                     |                              |
| $g=2\pi\hbar^2a/\mu$ | 赝势强度 | 能量·长度$^3$                                              |                          |                              |
| $\eta=\partial_r(r\Psi)|_{r\to0}$                                                               | 赝势的"正则部分"         | $1/$长度$^{3/2}$（约定） |
| $\psi_s(z)$                                        | 第$s$ 个横向模式（$n=2s$）的纵向波函数 | $1/$长度               |                              |
| $\kappa_s=(2/a_\perp)\sqrt{s+\epsilon}$                                                         | 闭道衰减常数             | 1/长度                       |
| $\epsilon=-k_z^2a_\perp^2/4$                                                                    | 无量纲能量参数           | —                           |
| $\Lambda[\xi,\epsilon]$                                                                         | 闭道求和函数，式 (32)    | —                           |
| $\mathcal C=-\zeta(1/2)=1.4603\ldots$                                                           | 由闭道求和产生的常数     | —                           |
| $\bar{\mathcal L}(\epsilon)$                                                                    | 式 (38) 的解析函数       | —                           |
| $f_{\rm even},f_{\rm odd}$                                                                      | 一维偶/奇散射振幅        | 长度                         |
| $a_{\rm 1D}$                                                                                    | 一维散射长度，式 (42)    | 长度                         |
| $g_{\rm 1D}=-\hbar^2/(\mu a_{\rm 1D})$                                                          | 有效一维耦合常数         | 能量·长度                   |
| $\Delta$                                                                                        | 一维偶道散射相位         | —                           |
| $\mathcal T$                                                                                    | 透射系数                 | —                           |
| $\Psi^{b},\Psi^{f}$                                                                             | Tonks 气体基态/费米基态  | —                           |
| $k_{\rm Fermi}$                                                                                 | 一维费米动量             | 1/长度                       |
| $\rho=N/L$                                                                                      | 一维线密度               | 1/长度                       |
| $N^*=L/(\pi|a_{\rm 1D}|)$                                                                       | 原子数上限               | —                           |
| $w(k_z)$                                                                                        | 一体动量分布             | 长度                         |

---

## 附录 E：本材料公式与原文公式对照

| 本材料     | 内容                                                       | 原文     |
| ---------- | ---------------------------------------------------------- | -------- |
| (1)(2)(3)  | 两体分离、相对哈密顿量、$a_\perp$                        | (2)(3)   |
| (4)(5)     | 赝势、$g=2\pi\hbar^2a/\mu$                               | (1)      |
| (6)        | 低能散射振幅与散射长度                                     | —       |
| (7)(8)     | 三维 Green 函数、Lippmann–Schwinger 解                    | —       |
| (9)–(13)  | 横向谐振子、$E_n$、$\phi_{0,0}$、$|\phi_{n,0}(0)|^2$ | (3)      |
| (14)–(17) | 能量、约束条件、渐近形式                                   | (4)(5)   |
| (18)–(30) | 模式展开、投影方程、Green 函数、$\psi_0,\psi_s$          | (6)(7)   |
| (31)(32)   | $\Psi(z,0)$、$\Lambda$ 定义                            | (8)      |
| (33)       | $\Lambda$ 展开、$\mathcal C,\bar{\mathcal L}$          | (9)(10)  |
| (34)–(42) | 自洽条件、$f_{\rm even}$、$a_{\rm 1D}$                 | (11)(12) |
| (43)       | $a_{\rm 1D}$ 的相位定义                                  | 原文文字 |
| (44)–(47) | 一维$\delta$ 势、$f^\delta$、透射系数                  | (13)     |
| (48)       | $g_{\rm 1D}$                                             | (14)     |
| (49)       | CIR 位置                                                   | —       |
| (50)       | 不可穿透条件                                               | (15)     |
| (51)–(53) | 多体$\delta$ 气体、硬核、Girardeau 映射                  | (16)     |
| (54)(55)   | Slater 行列式、$k_{\rm Fermi}$                           | (17)(18) |
| (56)       | 动量分布小$k$ 峰                                         | 图 2     |
| (57)       | 原子数上限                                                 | (19)     |
| (58)       | 适用条件$|a|\gg r_0$                                     | 原文文字 |

---

## 附录 F：自检清单（检验是否真正掌握）

1. 为什么质心分离后横向陷阱频率不变、但质量变成 $\mu=m/2$？$a_\perp$ 为什么是 $\sqrt{2\hbar/(m\omega_\perp)}$？
2. 赝势 (4) 中"正规化算符" $\frac{\partial}{\partial r}(r\,\cdot\,)$ 去掉了什么发散？$g=2\pi\hbar^2a/\mu$ 是怎么从 Lippmann–Schwinger 方程定出的？
3. 为什么只有 $m=0$ 的横向模式参与？为什么闭道求和从 $n=2$ 开始？
4. 为什么投影响应是一个严格的点源 $\eta\delta(z)$，而不是近似？
5. 验证 Green 函数 (27)(28) 的跳变条件。
6. $\Lambda$ 为什么在 $\xi\to0$ 不能逐项求导？$\mathcal C=-\zeta(1/2)$ 从哪一步冒出来？
7. 自洽条件 $\eta=\partial_z[z\Psi]|_{0+}$ 中，为什么恰是 $\Lambda$ 的 $2/\xi$ 项贡献了发散、$\mathcal L$ 项贡献了有限常数？
8. 从 (39)(40) 推导 (41)(42)（这是必须能亲手完成的代数）。
9. 一维 $\delta$ 势为什么零能全反射？$g_{\rm 1D}$ 与 $a_{\rm 1D}$ 的关系是什么？
10. CIR 在 $a=a_\perp/\mathcal C$ 处；为什么说它"不需要三维共振"？
11. Girardeau 映射为什么只在 1D 成立？$|\Psi^f|$ 满足哪三个条件使得它成为硬核玻色子的本征态？
12. $N^*\propto L\omega_\perp$ 的来源；代入 Rb/Na 的数字验证。
