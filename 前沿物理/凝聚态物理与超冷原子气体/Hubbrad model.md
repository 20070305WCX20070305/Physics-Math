# Hubbard Model and Mott Insulator

> **Context:** These notes assume a basic foundation in quantum mechanics, but not yet a full background in quantum many-body physics.

---

# 1. Why Do We Need the Hubbard Model?

A central problem in condensed matter physics is understanding a large number of interacting quantum particles.

For electrons in a solid, two effects are especially important:

1. **Electrons can move through the crystal.**
2. **Electrons interact with each other.**

The Hubbard model is a simplified model that captures the competition between these two effects.

The central question is:

> **What happens when quantum particles can hop between lattice sites while strongly interacting when they occupy the same site?**

The model is simple enough to write down but rich enough to describe important many-body phenomena such as:

* metallic behavior,
* Mott insulating behavior,
* strong correlations,
* quantum magnetism,
* quantum phase transitions,
* and potentially superconductivity.

---

# 2. From a Real Material to a Lattice Model

Consider electrons in a crystal.

The microscopic Hamiltonian can be written schematically as

$$
H =
\sum_i
\left[
\frac{\mathbf p_i^2}{2m}
+
V_{\rm crystal}(\mathbf r_i)
\right]
+
\sum_{i<j}
V_{\rm int}(\mathbf r_i-\mathbf r_j).
$$

Here:

* \(\frac{\mathbf p^2}{2m}\) is the kinetic energy,
* \(V_{\rm crystal}\) is the periodic potential produced by the crystal,
* \(V_{\rm int}\) describes electron-electron interactions.

This is a continuous-space description.

However, if the periodic potential is sufficiently strong, electrons can be approximately described as occupying localized orbitals around lattice sites.

This leads to a **lattice model**.

---

# 3. Optical Lattice and Ultracold Atoms

The same idea can be implemented experimentally using ultracold atoms.

Counter-propagating lasers can create a periodic potential such as

$$
V(x)=V_0\sin^2(kx).
$$

The potential contains many wells:

```text
Potential

 V
 │       /\        /\        /\        /\
 │      /  \      /  \      /  \      /  \
 │_____/    \____/    \____/    \____/    \____ x
       1         2         3         4
```

Atoms are localized around the minima of the potential.

Therefore, an ultracold-atom system can be engineered to behave like a lattice many-body system.

This is one reason ultracold atoms are useful for studying condensed-matter physics.

---

# 4. The Tight-Binding Picture

Let

$$
|w_i\rangle
$$

denote a localized orbital centered around lattice site \(i\). These are commonly called **Wannier states**.

The many-body wavefunction can be expanded in terms of these localized states.

The Hamiltonian contains matrix elements such as

$$
\langle w_i|H|w_j\rangle.
$$

For neighboring sites, we define approximately

$$
\boxed{
t=-\langle w_i|H|w_j\rangle
}
$$

where \(t\) is the **hopping amplitude**.

Therefore, the original kinetic energy

$$
\frac{p^2}{2m}
$$

does not literally disappear.

Instead, after projecting the continuous system onto localized lattice states, its effect is represented by the **hopping term**.

Schematically,

$$
\boxed{
\frac{p^2}{2m}+V_{\rm lattice}
\quad
\longrightarrow
\quad
\text{hopping with amplitude }t.
}
$$

This is an important connection between the microscopic system and the Hubbard model.

---

# 5. The Hubbard Hamiltonian

The standard Hubbard Hamiltonian is

$$
\boxed{
H=
-t
\sum_{\langle i,j\rangle,\sigma}
\left(
c^\dagger_{i\sigma}c_{j\sigma}
+
c^\dagger_{j\sigma}c_{i\sigma}
\right)
+
U\sum_i
n_{i\uparrow}n_{i\downarrow}
}
$$

or, using the common abbreviation,

$$
\boxed{
H=
-t
\sum_{\langle i,j\rangle,\sigma}
\left(
c^\dagger_{i\sigma}c_{j\sigma}
+h.c.
\right)
+
U\sum_i
n_{i\uparrow}n_{i\downarrow}.
}
$$

The model contains two essential energy scales:

$$
\boxed{t=\text{hopping energy}}
$$

and

$$
\boxed{U=\text{on-site interaction energy}}.
$$

The central dimensionless parameter is

$$
\boxed{\frac{U}{t}}.
$$

---

# 6. Meaning of the Operators

The Hubbard model is usually written in **second quantization**.

### Creation operator

$$
c^\dagger_{i\sigma}
$$

creates a particle at site \(i\) with spin \(\sigma\).

### Annihilation operator

$$
c_{i\sigma}
$$

removes a particle from site \(i\) with spin \(\sigma\).

### Number operator

$$
\boxed{
n_{i\sigma}
=
c^\dagger_{i\sigma}c_{i\sigma}
}
$$

counts particles with spin \(\sigma\) at site \(i\).

For spin-\(\frac12\) particles,

$$
\sigma=\uparrow,\downarrow.
$$

Therefore,

$$
n_{i\uparrow}n_{i\downarrow}
$$

is nonzero when site \(i\) is doubly occupied.

---

# 7. What Does "h.c." Mean?

The abbreviation

$$
\boxed{h.c.}
$$

means **Hermitian conjugate**.

For example,

$$
c^\dagger_i c_j+h.c.
$$

means

$$
c^\dagger_i c_j
+
c^\dagger_j c_i.
$$

Physically:

$$
c^\dagger_i c_j
$$

describes hopping

$$
j\rightarrow i,
$$

while

$$
c^\dagger_jc_i
$$

describes the reverse process

$$
i\rightarrow j.
$$

The Hamiltonian must be Hermitian so that its energy eigenvalues are real.

In general,

$$
A+h.c.
$$

means

$$
A+A^\dagger.
$$

---

# 8. Physical Meaning of the Two Terms

Write

$$
H=H_t+H_U.
$$

## 8.1 Hopping term

$$
H_t=
-t\sum_{\langle i,j\rangle,\sigma}
(c^\dagger_{i\sigma}c_{j\sigma}+h.c.).
$$

This describes particles moving between neighboring sites.

Conceptually,

```text
●     ○
 \___/
   t
```

Large \(t\) means strong tunneling and high mobility.

---

## 8.2 Interaction term

$$
H_U=
U\sum_i n_{i\uparrow}n_{i\downarrow}.
$$

This gives an energy cost \(U\) when two opposite-spin particles occupy the same site.

```text
Single occupation:       Double occupation:

    ↑                         ↑↓
    ●                         ●
```

For repulsive interactions,

$$
U>0.
$$

Large \(U\) strongly suppresses double occupation.

---

# 9. The Central Competition: \(U\) vs. \(t\)

The Hubbard model is fundamentally about the competition

$$
\boxed{
\text{kinetic energy}
\quad\leftrightarrow\quad
\text{interaction energy}
}
$$

or

$$
\boxed{
t
\quad\leftrightarrow\quad
U.
}
$$

## Weakly interacting regime

$$
U\ll t.
$$

Particles can move relatively easily.

The system is more weakly correlated.

---

## Strongly interacting regime

$$
U\gg t.
$$

Particles strongly avoid double occupation.

Their motion becomes strongly correlated.

This is the regime where Mott physics becomes important.

---

# 10. What Are Strong Correlations?

In a weakly interacting system, we may approximately think about particles individually.

For example:

```text
particle 1 → moves
particle 2 → moves
particle 3 → moves
```

But when

$$
U\gg t,
$$

the particles strongly affect each other.

The motion of particle 1 depends on whether neighboring sites are occupied.

Thus,

$$
\boxed{
\text{the particles cannot be treated independently}.
}
$$

This is the essence of **strongly correlated many-body physics**.

---

# 11. Mott Insulator

A **Mott insulator** is an insulating state produced primarily by strong particle-particle interactions.

Consider a lattice with approximately one particle per site:

```text
●     ●     ●     ●     ●
```

Suppose

$$
U\gg t.
$$

For a particle to hop to a neighboring occupied site,

```text
Before:

●     ●

After:

○     ●●
```

the intermediate state has double occupation.

That costs energy

$$
U.
$$

Therefore, real charge motion is strongly suppressed.

The result is approximately

$$
\boxed{
\text{one particle per site}
+
U\gg t
\quad\Rightarrow\quad
\text{Mott insulating behavior}.
}
$$

---

# 12. Why Is a Mott Insulator Different from an Ordinary Insulator?

An ordinary band insulator can be understood largely from single-particle band structure.

Very schematically,

$$
\text{filled band}
+
\text{energy gap}
\rightarrow
\text{insulator}.
$$

A Mott insulator is fundamentally different.

Its insulating behavior arises from interactions:

$$
\boxed{
\text{strong interaction}
\rightarrow
\text{suppressed charge motion}
\rightarrow
\text{insulator}.
}
$$

Thus, a Mott insulator is an example of an **emergent many-body phenomenon**.

---

# 13. Are the Particles Completely Frozen?

No.

This is an important subtlety.

In a Mott insulator,

$$
t\neq0.
$$

The particles can still undergo **virtual hopping**.

For example,

```text
Initial:

↑       ↓

        ↓

Virtual state:

↑↓      ○

        ↓

Final:

↑       ↓
```

The intermediate state costs energy \(U\), so the process is suppressed but not completely absent.

This virtual hopping becomes extremely important for magnetism.

---

# 14. From the Hubbard Model to the Heisenberg Model

Suppose each site contains one spin-\(\frac12\) particle.

Each particle can be in

$$
|\uparrow\rangle
$$

or

$$
|\downarrow\rangle.
$$

Although real charge motion is suppressed, virtual hopping allows neighboring spins to interact.

For

$$
U\gg t,
$$

the effective low-energy Hamiltonian becomes approximately

$$
\boxed{
H_{\rm eff}
=
J
\sum_{\langle i,j\rangle}
\mathbf S_i\cdot\mathbf S_j
}
$$

with

$$
\boxed{
J\approx\frac{4t^2}{U}.
}
$$

This is the **Heisenberg model**.

Therefore,

$$
\boxed{
\text{Hubbard model}
\xrightarrow{U\gg t}
\text{Heisenberg spin model}.
}
$$

The important physical process is

$$
\boxed{
\text{virtual hopping}
\rightarrow
\text{effective spin-spin interaction}.
}
$$

This is the origin of quantum magnetism in this simple picture.

---

# 15. Charge vs. Spin: An Important Distinction

A Mott insulator can have very little **charge mobility** while still having significant **spin dynamics**.

### Charge sector

Real hopping is energetically expensive:

$$
\boxed{
\text{charge motion}\quad\text{suppressed}
}
$$

leading to insulating behavior.

### Spin sector

Virtual hopping produces

$$
J\sim\frac{4t^2}{U},
$$

so spins can still interact:

$$
\boxed{
\text{spin dynamics}\quad\text{remain possible}.
}
$$

Thus:

$$
\boxed{
\text{Mott insulator}
=
\text{charge localization}
+
\text{potentially rich spin physics}.
}
$$

---

# 16. Why Ultracold Atoms Are Useful

Ultracold atoms provide a highly controllable realization of the Hubbard model.

A typical experimental route is

$$
\boxed{
\text{ultracold atoms}
\rightarrow
\text{optical lattice}
\rightarrow
\text{localized atoms}
\rightarrow
\text{Hubbard model}.
}
$$

The experimental parameters can be controlled to change the effective model parameters.

For example:

### Optical lattice depth

Changing the lattice depth changes the tunneling amplitude:

$$
V_0\uparrow
\quad\Rightarrow\quad
t\downarrow.
$$

### Atomic interaction

The interaction strength can be tuned experimentally, changing

$$
U.
$$

Therefore, researchers can explore different values of

$$
\boxed{\frac{U}{t}}.
$$

This allows experimental investigation of strongly correlated regimes.

---

# 17. The Conceptual Map

The whole story can be summarized as:

```text
                 Ultracold atoms
                       │
                       ▼
                Optical lattice
                       │
                       ▼
             Localized Wannier states
                       │
                       ▼
                Hubbard model
                       │
             ┌─────────┴─────────┐
             ▼                   ▼
          hopping              interaction
             t                     U
             │                   │
             └─────────┬─────────┘
                       ▼
                      U/t
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
       weak U       strong U      intermediate
          │            │
          ▼            ▼
      mobile       Mott regime
     particles         │
                       ▼
              suppressed charge
                       │
                       ▼
               Mott insulator
                       │
                       ▼
                virtual hopping
                       │
                       ▼
              spin-spin coupling
                       │
                       ▼
              Heisenberg model
                       │
                       ▼
              Quantum magnetism
```

---

# 18. Key Equations to Remember

### Microscopic Hamiltonian

$$
H=
\frac{p^2}{2m}
+
V_{\rm lattice}
+
V_{\rm int}.
$$

### Optical lattice

$$
V(x)=V_0\sin^2(kx).
$$

### Hubbard Hamiltonian

$$
\boxed{
H=
-t\sum_{\langle i,j\rangle,\sigma}
(c^\dagger_{i\sigma}c_{j\sigma}+h.c.)
+
U\sum_i n_{i\uparrow}n_{i\downarrow}.
}
$$

### Number operator

$$
\boxed{
n_{i\sigma}=c^\dagger_{i\sigma}c_{i\sigma}.
}
$$

### Important energy ratio

$$
\boxed{
\frac{U}{t}.
}
$$

### Strongly correlated regime

$$
\boxed{
U\gg t.
}
$$

### Effective magnetic interaction

$$
\boxed{
J\approx\frac{4t^2}{U}.
}
$$

---

# 19. The Most Important Physical Picture

If you remember only one picture, remember this:

$$
\boxed{
\text{Hubbard model}
=
\underbrace{\text{particles want to move}}_{t}
+
\underbrace{\text{particles interact}}_{U}
}
$$

The competition between these two tendencies produces rich many-body physics.

In particular,

$$
\boxed{
U\gg t
\Rightarrow
\text{strong correlations}
\Rightarrow
\text{suppressed charge motion}
\Rightarrow
\text{Mott insulator}.
}
$$

But because \(t\neq0\),

$$
\boxed{
\text{virtual hopping}
\Rightarrow
J\sim\frac{4t^2}{U}
\Rightarrow
\text{spin interactions}.
}
$$

Therefore, a Mott insulator is not simply a collection of completely frozen particles. It is a **strongly correlated quantum state with localized charge but potentially rich spin dynamics**.

---

# 20. A Useful Learning Roadmap

A natural progression for studying this subject is:

$$
\boxed{
\text{Quantum mechanics}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Identical particles + spin}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Second quantization}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Tight-binding approximation}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Hubbard model}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Mott insulator}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Heisenberg model + quantum magnetism}
}
$$

$$
\downarrow
$$

$$
\boxed{
\text{Quantum phase transitions / strongly correlated systems}
}
$$

This sequence provides a very natural bridge from undergraduate quantum mechanics to **ultracold atoms and modern quantum many-body physics**.
