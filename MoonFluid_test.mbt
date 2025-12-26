// MoonFluid - 测试代码

///|
test "hydrostatic pressure" {
  // 水深 10m，水密度 1000 kg/m³
  // p = p0 + ρgh = 101325 + 1000*9.80665*10 = 199391.5 Pa
  let p = hydrostatic_pressure_std(ATMOSPHERIC_PRESSURE, RHO_WATER, 10.0)
  assert_eq(approx_equal(p, 199391.5, 50.0), true)

  // 浮力
  let fb = buoyant_force_std(RHO_WATER, 1.0) // 1 m³ 水
  assert_eq(approx_equal(fb, 9806.65, 1.0), true)
}

///|
test "continuity equation" {
  // 体积流量
  let q = volume_flow_rate(0.1, 2.0) // 0.1 m² * 2 m/s
  assert_eq(q, 0.2)

  // 流速计算
  let v = velocity_from_flow_rate(0.2, 0.1)
  assert_eq(v, 2.0)

  // 质量流量
  let mdot = mass_flow_rate(RHO_WATER, 0.2)
  assert_eq(mdot, 200.0)
}

///|
test "bernoulli equation" {
  // 总压头: z=5, p/(ρg)=10.19, v²/(2g)=0.459
  let h = bernoulli_total_head_std(5.0, 100000.0, RHO_WATER, 3.0)
  assert_eq(approx_equal(h, 15.65, 0.5), true)

  // 速度水头
  let hv = velocity_head(5.0, G)
  assert_eq(approx_equal(hv, 1.274, 0.01), true)
}

///|
test "reynolds number" {
  // Re = ρvD/μ
  // 水: ρ=1000, v=1 m/s, D=0.1 m, μ=0.001002
  // Re = 1000*1*0.1/0.001002 ≈ 99800
  let re = reynolds_number(RHO_WATER, 1.0, 0.1, MU_WATER_20C)
  assert_eq(approx_equal(re, 99800.0, 500.0), true)

  // 判断流动状态
  assert_eq(is_turbulent_pipe(5000.0), true)
  assert_eq(is_turbulent_pipe(1000.0), false)
}

///|
test "laminar flow" {
  // 层流流量: Q = πΔpR⁴/(8μL)
  // Δp=1000 Pa, μ=0.001, L=10 m, R=0.05 m
  // Q = π*1000*0.05⁴/(8*0.001*10) = π*1000*0.00000625/0.08 = 0.245 m³/s
  let q = laminar_flow_rate(1000.0, MU_WATER_20C, 10.0, 0.05)
  assert_eq(approx_equal(q, 0.245, 0.01), true)

  // 层流摩擦系数
  let f = laminar_friction_factor(1000.0)
  assert_eq(f, 0.064)
}

///|
test "darcy weisbach" {
  // 摩擦损失水头
  // f=0.02, L=100 m, D=0.2 m, v=2 m/s
  // hf = 0.02*(100/0.2)*(4/(2*9.81)) = 0.02*500*0.204 = 2.04 m
  let hf = darcy_weisbach_head_loss_std(0.02, 100.0, 0.2, 2.0)
  assert_eq(approx_equal(hf, 2.04, 0.1), true)
}

///|
test "local losses" {
  // 局部损失
  let hl = local_head_loss_std(0.5, 3.0)
  assert_eq(approx_equal(hl, 0.229, 0.01), true)

  // 突然扩大
  let k_exp = sudden_expansion_loss_coefficient(0.1, 0.2)
  assert_eq(approx_equal(k_exp, 0.25, 0.01), true)
}

///|
test "manning formula" {
  // 矩形明渠: b=2m, h=1m, n=0.013, S=0.001
  // A = 2*1 = 2 m², P = 2+2*1 = 4 m, R = 2/4 = 0.5 m
  // Q = (1/0.013)*2*0.5^(2/3)*0.001^(1/2)
  let r = hydraulic_radius_rect(2.0, 1.0)
  assert_eq(approx_equal(r, 0.5, 0.01), true)
  let q = manning_flow_rate(0.013, 2.0, 0.5, 0.001)
  assert_eq(q > 0.0, true)
}

///|
test "froude number" {
  // Fr = v/sqrt(gh)
  // v=2 m/s, h=1 m
  // Fr = 2/sqrt(9.81*1) = 0.639
  let fr = froude_number(2.0, G, 1.0)
  assert_eq(approx_equal(fr, 0.639, 0.01), true)

  // 流动状态
  assert_eq(flow_regime(0.5), "subcritical")
  assert_eq(flow_regime(1.0), "critical")
  assert_eq(flow_regime(2.0), "supercritical")
}

///|
test "drag and lift" {
  // 阻力: F_D = 0.5*ρ*v²*A*C_D
  // ρ=1.204, v=10 m/s, A=1 m², C_D=0.5
  // F_D = 0.5*1.204*100*1*0.5 = 30.1 N
  let fd = drag_force(RHO_AIR, 10.0, 1.0, 0.5)
  assert_eq(approx_equal(fd, 30.1, 0.1), true)

  // 升力: F_L = 0.5*ρ*v²*A*C_L = 0.5*1.204*100*1*0.8 = 48.16 N
  let fl = lift_force(RHO_AIR, 10.0, 1.0, 0.8)
  assert_eq(approx_equal(fl, 48.16, 0.1), true)
}

///|
test "dimensionless numbers" {
  // 韦伯数
  let we = weber_number(RHO_WATER, 5.0, 0.1, SIGMA_WATER)
  assert_eq(we > 0.0, true)

  // 马赫数
  let ma = mach_number(340.0, 340.0)
  assert_eq(ma, 1.0)

  // 欧拉数: Eu = Δp/(ρv²) = 1000/(1000*25) = 0.04
  let eu = euler_number(1000.0, RHO_WATER, 5.0)
  assert_eq(approx_equal(eu, 0.04, 0.01), true)
}

///|
test "tools functions" {
  // 平方、立方
  assert_eq(square(3.0), 9.0)
  assert_eq(cube(2.0), 8.0)

  // 立方根
  assert_eq(approx_equal(cbrt(27.0), 3.0, 0.0001), true)
  assert_eq(approx_equal(cbrt(8.0), 2.0, 0.0001), true)

  // 对数 (近似值，容差较大)
  assert_eq(approx_equal(log10_approx(100.0), 2.0, 0.3), true)
}
