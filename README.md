# Polyhedral Omega for Optimization

> Graduation Project – Department of Computer Engineering, Gebze Technical University  
> **Author:** Kevser Yolcu  
> **Supervisor:** Asst. Prof. Zafeirakis Zafeirakopoulos  
> **Date:** January 2021  
> **Location:** Gebze, KOCAELİ

## 📌 Abstract

Polyhedral Omega is a method for solving Linear Diophantine Systems (LDS) of inequalities. It computes all solutions to the system and can be adapted into an Integer Linear Programming (ILP) method using binary search strategies.

## 📚 Table of Contents

- [Abstract](#-abstract)
- [Introduction](#-introduction)
- [Method](#-method)
- [Results](#-results)
- [Discussion](#-discussion)
- [References](#-references)

## 🔍 Introduction

### Linear Diophantine Systems and Rational Functions

A **Linear Diophantine Equation** is an equation of the form:

```
ax + by = c, where a, b, c, x, y ∈ ℤ
```

A **Linear Diophantine System (LDS)** generalizes this into systems:

```
Ax = b, where A is an integer matrix, x is a vector of unknowns
```

A **rational function** is defined as:

```
f(x) = P(x)/Q(x), where P and Q are polynomials and Q ≠ 0
```

### Symbolic Cones

A symbolic cone is defined as a triple:

```
C = (V, q, o)
```

Where:
- `V ∈ ℤ^(d+m)×d` is the generator matrix
- `q ∈ ℤ^(d+m)` is the apex vector
- `o ∈ {0,1}^d` is the openness vector

### Julia Language

Julia is a high-performance, dynamic programming language ideal for numerical computing.

## 🧠 Method

### 1. MacMahon Lifting

Transforms input constraints into a symbolic cone by:
- Extending identity matrix `I_d` with matrix `A`
- Creating apex vector from zero vector and `-b`
- Generating openness vector

### 2. Iterative Elimination

Eliminates variables one-by-one by:
- Intersecting the cone with half-spaces
- Projecting out last coordinates
- Producing signed cones for further processing

### 3. Conversion to Rational Function

Uses **Smith Normal Form** to enumerate lattice points in the fundamental parallelepiped:

```
V = USW
```

Where:
- `U`, `W` are unimodular matrices
- `S` is a diagonal matrix

### 4. Normal Form of Rational Function

If solution space is finite, the output is a polynomial listing all integer solutions.

## 📈 Results

**Input:**

```text
A = [ [1, -1],
      [-2, 1] ],
b = [0, 0]
```

**After MacMahon Lifting:**

```text
V = [
  [1, 0],
  [0, 1],
  [1, -1],
  [-2, 1]
]
q = [0, 0, 0, 0]
o = [false, false]
```

**After Iterative Elimination:**

```text
Cone 1:
V1 = [[1, 0], [2, 1]]
q1 = [0, 0]
o1 = [true, true]
sgn1 = +1

Cone 2:
V2 = [[1, 1], [2, 1]]
q2 = [0, 0]
o2 = [true, false]
sgn2 = +1

Cone 3:
V3 = [[1, 0], [1, 1]]
q3 = [0, 0]
o3 = [true, true]
sgn3 = -1
```

**Lattice Points (Solutions):**

```text
p1 = [1, 3]
p2 = [1, 2]
p3 = [1, 2]
```

## 💬 Discussion

The project successfully implemented the **Polyhedral Omega** algorithm in Julia, demonstrating:
- Construction of symbolic cones from constraints
- Iterative variable elimination
- Rational function generation via Smith Normal Form

Although functional, further optimization of the algorithm would require a longer development cycle.

## 📚 References

1. [Wikipedia - Rational Function](https://en.wikipedia.org/wiki/Rational_function)  
2. Martin Heller, *What is Julia?*, InfoWorld, 2018  
3. Breuer & Zafeirakopoulos, *Polyhedral Omega: A New Algorithm for Solving Linear Diophantine Systems*
