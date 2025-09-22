# Simulation and Modeling of Dynamic Systems Projects

This repository contains the implementation of three projects in **Simulation and Modeling of Dynamic Systems**, focusing on parameter estimation, real-time adaptive methods, and robust modeling approaches.  
The work was developed as part of the *Simulation and Modeling of Dynamic Systems* course at the **Department of Electrical and Computer Engineering, Aristotle University of Thessaloniki.**

--- 

## 📌 Projects Overview

### **Project 1 – Parameter Estimation with Least Squares**
- System: Simple pendulum with torque input (linearized dynamics).  
- Tasks:
  - Derivation of state-space equations and transfer function.
  - Simulation of the system response using ODE solvers in MATLAB.
  - Parameter estimation (`m`, `L`, `c`) via the **Least Squares method**:
    - Case 1: full state vector measurable.
    - Case 2: only the angle `q(t)` measurable.
  - Study of robustness under:
    - measurement noise,
    - different sampling periods `Ts`,
    - varying input amplitude `A0`.

---

### **Project 2 – Real-Time Parameter Estimation Methods**
- Systems:  
  1. **Mass-spring-damper** with external force.  
  2. **Roll angle dynamics** of an aircraft (nonlinear system).  
- Tasks:
  - Real-time estimators using the **Gradient method** and **Lyapunov-based methods** (parallel and mixed structures).
  - Performance evaluation under:
    - constant and sinusoidal inputs,
    - measurement noise,
    - external disturbances.
  - Comparative analysis of parameter convergence and system tracking accuracy.

---

### **Final Project – Robust Real-Time Methods and Model Selection**
- Part 1: **Linear system with unknown matrices A, B**  
  - Real-time estimation of system matrices.  
  - Robust extensions with bias/disturbance modeling.  
  - Analysis of stability and accuracy under varying bias levels.  
- Part 2: **Nonlinear system approximation**  
  - Model structure selection using basis functions (e.g., polynomial, Gaussian).  
  - Real-time parameter estimation and cross-validation for model comparison.  
  - Final model selection based on trade-off between accuracy and complexity.  
  - Stability evaluation with independent test datasets.
