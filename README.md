# MoonFluid - Fluid Mechanics Calculation Library

**MoonFluid** is a comprehensive fluid mechanics calculation library written in **Moonbit**, providing functions for hydrostatics, fluid dynamics, pipe flow, open channel flow, and aerodynamic calculations.

---

## Features

### Fluid Statics

* **Hydrostatic Pressure**: Pressure calculation in static fluids (p = p₀ + ρgh)
* **Buoyant Force**: Buoyancy calculations (F_b = ρgV)
* **Total Pressure**: Pressure on submerged surfaces
* **Pressure Center**: Center of pressure on submerged planes
* **Manometry**: Fluid pressure measurement calculations

---

### Continuity and Bernoulli

* **Continuity Equation**: Mass and volume flow rate conservation (Q = A₁v₁ = A₂v₂)
* **Bernoulli Equation**: Energy conservation in fluid flow (p/ρg + v²/2g + z = constant)
* **Energy Head**: Total head calculations in flow systems
* **Velocity Head**: Kinetic energy head calculations

---

### Flow Characteristics

* **Reynolds Number**: Flow regime classification (Re = ρvD/μ)
* **Flow Regime Analysis**: Laminar vs turbulent flow determination
* **Froude Number**: Open channel flow regime analysis (Fr = v/√(gh))
* **Critical Flow**: Critical flow conditions and analysis

---

### Pipe Flow

* **Laminar Flow**: Hagen-Poiseuille flow calculations
* **Turbulent Flow**: Darcy-Weisbach friction loss calculations
* **Local Losses**: Sudden expansion, contraction, and fitting losses
* **Total Head Loss**: Combined friction and local losses
* **Pipe Sizing**: Flow rate and diameter calculations

---

### Open Channel Flow

* **Manning Equation**: Open channel flow calculations
* **Hydraulic Radius**: Channel geometry calculations
* **Normal Depth**: Uniform flow depth calculations
* **Critical Depth**: Critical flow depth analysis

---

### Aerodynamics and Forces

* **Drag Force**: Fluid resistance calculations (F_D = ½ρv²AC_D)
* **Lift Force**: Aerodynamic lift calculations (F_L = ½ρv²AC_L)
* **Drag Coefficient**: Flow resistance coefficients
* **Lift Coefficient**: Aerodynamic lift coefficients

---

### Dimensionless Numbers

* **Weber Number**: Surface tension effects (We = ρv²L/σ)
* **Mach Number**: Compressible flow analysis (Ma = v/a)
* **Euler Number**: Pressure difference analysis (Eu = Δp/ρv²)
* **Dimensionless Analysis**: Flow similarity and scaling

---

### Mathematical Utilities

* **Power Functions**: Square, cube, cube root calculations
* **Logarithmic Functions**: Base-10 logarithm approximations
* **Trigonometric Functions**: Advanced angle calculations
* **Numerical Methods**: Approximation and interpolation functions

---

## Usage

### Fluid Statics

```moonbit
test "hydrostatic pressure" {
  // Water depth 10m, water density 1000 kg/m³
  // p = p0 + ρgh = 101325 + 1000*9.80665*10 = 199391.5 Pa
  let p = hydrostatic_pressure_std(ATMOSPHERIC_PRESSURE, RHO_WATER, 10.0)
  assert_eq(approx_equal(p, 199391.5, 50.0), true)

  // Buoyant force
  let fb = buoyant_force_std(RHO_WATER, 1.0) // 1 m³ water
  assert_eq(approx_equal(fb, 9806.65, 1.0), true)
}
```

---

### Continuity Equation

```moonbit
test "continuity equation" {
  // Volume flow rate
  let q = volume_flow_rate(0.1, 2.0) // 0.1 m² * 2 m/s
  assert_eq(q, 0.2)

  // Velocity calculation from flow rate
  let v = velocity_from_flow_rate(0.2, 0.1)
  assert_eq(v, 2.0)

  // Mass flow rate
  let mdot = mass_flow_rate(RHO_WATER, 0.2)
  assert_eq(mdot, 200.0)
}
```

---

### Bernoulli Equation

```moonbit
test "bernoulli equation" {
  // Total head: z=5, p/(ρg)=10.19, v²/(2g)=0.459
  let h = bernoulli_total_head_std(5.0, 100000.0, RHO_WATER, 3.0)
  assert_eq(approx_equal(h, 15.65, 0.5), true)

  // Velocity head
  let hv = velocity_head(5.0, G)
  assert_eq(approx_equal(hv, 1.274, 0.01), true)
}
```

---

### Reynolds Number

```moonbit
test "reynolds number" {
  // Re = ρvD/μ
  // Water: ρ=1000, v=1 m/s, D=0.1 m, μ=0.001002
  // Re = 1000*1*0.1/0.001002 ≈ 99800
  let re = reynolds_number(RHO_WATER, 1.0, 0.1, MU_WATER_20C)
  assert_eq(approx_equal(re, 99800.0, 500.0), true)

  // Determine flow regime
  assert_eq(is_turbulent_pipe(5000.0), true)
  assert_eq(is_turbulent_pipe(1000.0), false)
}
```

---

### Laminar Flow

```moonbit
test "laminar flow" {
  // Laminar flow rate: Q = πΔpR⁴/(8μL)
  // Δp=1000 Pa, μ=0.001, L=10 m, R=0.05 m
  // Q = π*1000*0.05⁴/(8*0.001*10) = π*1000*0.00000625/0.08 = 0.245 m³/s
  let q = laminar_flow_rate(1000.0, MU_WATER_20C, 10.0, 0.05)
  assert_eq(approx_equal(q, 0.245, 0.01), true)

  // Laminar friction factor
  let f = laminar_friction_factor(1000.0)
  assert_eq(f, 0.064)
}
```

---

### Darcy-Weisbach Equation

```moonbit
test "darcy weisbach" {
  // Friction head loss
  // f=0.02, L=100 m, D=0.2 m, v=2 m/s
  // hf = 0.02*(100/0.2)*(4/(2*9.81)) = 0.02*500*0.204 = 2.04 m
  let hf = darcy_weisbach_head_loss_std(0.02, 100.0, 0.2, 2.0)
  assert_eq(approx_equal(hf, 2.04, 0.1), true)
}
```

---

### Local Losses

```moonbit
test "local losses" {
  // Local head loss
  let hl = local_head_loss_std(0.5, 3.0)
  assert_eq(approx_equal(hl, 0.229, 0.01), true)

  // Sudden expansion loss coefficient
  let k_exp = sudden_expansion_loss_coefficient(0.1, 0.2)
  assert_eq(approx_equal(k_exp, 0.25, 0.01), true)
}
```

---

### Manning Formula

```moonbit
test "manning formula" {
  // Rectangular open channel: b=2m, h=1m, n=0.013, S=0.001
  // A = 2*1 = 2 m², P = 2+2*1 = 4 m, R = 2/4 = 0.5 m
  // Q = (1/0.013)*2*0.5^(2/3)*0.001^(1/2)
  let r = hydraulic_radius_rect(2.0, 1.0)
  assert_eq(approx_equal(r, 0.5, 0.01), true)
  let q = manning_flow_rate(0.013, 2.0, 0.5, 0.001)
  assert_eq(q > 0.0, true)
}
```

---

### Froude Number

```moonbit
test "froude number" {
  // Fr = v/sqrt(gh)
  // v=2 m/s, h=1 m
  // Fr = 2/sqrt(9.81*1) = 0.639
  let fr = froude_number(2.0, G, 1.0)
  assert_eq(approx_equal(fr, 0.639, 0.01), true)

  // Flow regime classification
  assert_eq(flow_regime(0.5), "subcritical")
  assert_eq(flow_regime(1.0), "critical")
  assert_eq(flow_regime(2.0), "supercritical")
}
```

---

### Drag and Lift Forces

```moonbit
test "drag and lift" {
  // Drag force: F_D = 0.5*ρ*v²*A*C_D
  // ρ=1.204, v=10 m/s, A=1 m², C_D=0.5
  // F_D = 0.5*1.204*100*1*0.5 = 30.1 N
  let fd = drag_force(RHO_AIR, 10.0, 1.0, 0.5)
  assert_eq(approx_equal(fd, 30.1, 0.1), true)

  // Lift force: F_L = 0.5*ρ*v²*A*C_L = 0.5*1.204*100*1*0.8 = 48.16 N
  let fl = lift_force(RHO_AIR, 10.0, 1.0, 0.8)
  assert_eq(approx_equal(fl, 48.16, 0.1), true)
}
```

---

### Dimensionless Numbers

```moonbit
test "dimensionless numbers" {
  // Weber number
  let we = weber_number(RHO_WATER, 5.0, 0.1, SIGMA_WATER)
  assert_eq(we > 0.0, true)

  // Mach number
  let ma = mach_number(340.0, 340.0)
  assert_eq(ma, 1.0)

  // Euler number: Eu = Δp/(ρv²) = 1000/(1000*25) = 0.04
  let eu = euler_number(1000.0, RHO_WATER, 5.0)
  assert_eq(approx_equal(eu, 0.04, 0.01), true)
}
```

---

### Mathematical Tools

```moonbit
test "tools functions" {
  // Square and cube
  assert_eq(square(3.0), 9.0)
  assert_eq(cube(2.0), 8.0)

  // Cube root
  assert_eq(approx_equal(cbrt(27.0), 3.0, 0.0001), true)
  assert_eq(approx_equal(cbrt(8.0), 2.0, 0.0001), true)

  // Logarithm (approximate, with larger tolerance)
  assert_eq(approx_equal(log10_approx(100.0), 2.0, 0.3), true)
}
```

---

## Parameter Ranges

### Valid Input Ranges

* **Pressure**: 0–10⁸ Pa (fluid pressures)
* **Density**: 0.1–2000 kg/m³ (fluid densities)
* **Viscosity**: 10⁻⁶–10⁻² Pa·s (dynamic viscosities)
* **Velocity**: 0–100 m/s (flow velocities)
* **Length/Diameter**: 10⁻⁶–10² m (geometrical dimensions)
* **Flow Rate**: 10⁻⁹–10³ m³/s (volumetric flow rates)
* **Reynolds Number**: 1–10⁷ (flow regime classification)
* **Froude Number**: 0.1–10 (open channel flow)

---

### Typical Fluid Properties

* **Water (20°C)**: ρ = 998.2 kg/m³, μ = 1.002×10⁻³ Pa·s
* **Air (20°C)**: ρ = 1.204 kg/m³, μ = 1.81×10⁻⁵ Pa·s
* **Standard Gravity**: g = 9.80665 m/s²
* **Atmospheric Pressure**: p₀ = 101,325 Pa
* **Water Surface Tension**: σ = 0.0728 N/m
* **Speed of Sound in Air**: a = 343 m/s

---

## Testing

The project includes a comprehensive test suite covering all major functionalities:

```bash
moon test
```

### Test Coverage

* Fluid statics (pressure, buoyancy, pressure centers)
* Continuity and Bernoulli equations
* Flow regime analysis (Reynolds, Froude numbers)
* Pipe flow calculations (laminar, turbulent, losses)
* Open channel flow (Manning equation, hydraulic radius)
* Aerodynamic forces (drag, lift coefficients)
* Dimensionless numbers and similarity
* Mathematical utility functions
* Boundary conditions and edge cases

---

## Technical Details

### Physical Principles

* **Conservation Laws**: Mass, energy, and momentum conservation
* **Newtonian Fluids**: Viscosity and shear stress relationships
* **Pressure Fields**: Hydrostatic and hydrodynamic pressure distributions
* **Flow Similarity**: Dimensional analysis and scaling laws
* **Boundary Layers**: Flow behavior near solid surfaces

---

### Mathematical Methods

* **Algebraic Equations**: Direct calculations for statics and simple flows
* **Differential Equations**: Energy and momentum balance equations
* **Empirical Correlations**: Friction factors and loss coefficients
* **Numerical Approximations**: Root finding and interpolation methods
* **Statistical Analysis**: Uncertainty and error propagation

---

### Engineering Standards

* **SI Units**: Consistent use of International System of Units
* **Fluid Mechanics Conventions**: Standard notation and terminology
* **Engineering Accuracy**: Appropriate precision for engineering calculations
* **Safety Margins**: Conservative design practices

---

### Computational Considerations

* **Numerical Stability**: Robust algorithms for wide parameter ranges
* **Physical Limits**: Validation of input parameters against physical bounds
* **Unit Consistency**: Automatic dimensional analysis
* **Performance Optimization**: Efficient calculations for engineering applications

---

## Notes

1. **Units**: All calculations use SI units (meters, kilograms, seconds, Pascals)
2. **Fluid Properties**: Default values provided for common fluids (water, air)
3. **Flow Assumptions**: Incompressible, Newtonian fluid assumptions unless specified
4. **Temperature Effects**: Standard conditions assumed unless temperature specified
5. **Accuracy**: Engineering precision suitable for design calculations
6. **Boundary Conditions**: Physical limits enforced to prevent unrealistic results
7. **Convergence**: Iterative methods include convergence checks
8. **Validation**: Results should be verified against experimental data or established correlations

---

## Application Scenarios

* **Hydraulic Engineering**: Pipeline design, pump selection, pressure vessel analysis
* **Environmental Engineering**: Water distribution systems, wastewater treatment
* **Mechanical Engineering**: Fluid machinery design, heat exchangers, ventilation systems
* **Civil Engineering**: Open channel flow, drainage systems, hydraulic structures
* **Aerospace Engineering**: Aerodynamic analysis, wind tunnel testing, aircraft design
* **Chemical Engineering**: Process piping, mixing systems, fluid transport
* **Marine Engineering**: Ship hydrodynamics, offshore structures, underwater systems
* **Geotechnical Engineering**: Groundwater flow, seepage analysis, dam design
* **Energy Systems**: Turbomachinery, renewable energy, power plant design
* **Research Applications**: Fluid dynamics research, computational fluid dynamics validation

---

## Version Information

The current version (0.1.0) implements **comprehensive fluid mechanics calculation functions** including:

* Fluid statics (pressure, buoyancy, surface forces)
* Flow analysis (continuity, Bernoulli, energy equations)
* Flow characteristics (Reynolds, Froude, Mach numbers)
* Pipe flow (laminar, turbulent, friction losses)
* Open channel flow (Manning equation, critical flow)
* Aerodynamic forces (drag, lift, coefficients)
* Dimensionless analysis (Weber, Euler numbers)
* Mathematical utilities (power functions, logarithms)

The library is actively developed and aims to provide comprehensive coverage of fluid mechanics while being optimized for MoonBit.
