// MoonFluid - 流体力学库
// 主模块：包含核心计算功能

// ==================== 流体静力学 ====================

///|
/// 计算静水压强
/// p = p₀ + ρgh
/// 参数: p0 - 表面压强 (Pa), rho - 密度 (kg/m³), g - 重力加速度 (m/s²), h - 深度 (m)
/// 返回: 压强 (Pa)
pub fn hydrostatic_pressure(
  p0 : Double,
  rho : Double,
  g : Double,
  h : Double,
) -> Double {
  p0 + rho * g * h
}

///|
/// 计算静水压强 (使用标准重力)
pub fn hydrostatic_pressure_std(
  p0 : Double,
  rho : Double,
  h : Double,
) -> Double {
  hydrostatic_pressure(p0, rho, G, h)
}

///|
/// 计算浮力
/// F_b = ρgV
/// 参数: rho - 流体密度 (kg/m³), g - 重力加速度 (m/s²), volume - 排开体积 (m³)
/// 返回: 浮力 (N)
pub fn buoyant_force(rho : Double, g : Double, volume : Double) -> Double {
  rho * g * volume
}

///|
/// 计算浮力 (使用标准重力)
pub fn buoyant_force_std(rho : Double, volume : Double) -> Double {
  buoyant_force(rho, G, volume)
}

///|
/// 计算作用在平面上的总压力
/// F = ρgh_c * A
/// 参数: rho - 密度 (kg/m³), g - 重力加速度 (m/s²), h_c - 形心深度 (m), area - 面积 (m²)
/// 返回: 总压力 (N)
pub fn total_pressure_plane(
  rho : Double,
  g : Double,
  h_c : Double,
  area : Double,
) -> Double {
  rho * g * h_c * area
}

///|
/// 计算作用在平面上的总压力 (使用标准重力)
pub fn total_pressure_plane_std(
  rho : Double,
  h_c : Double,
  area : Double,
) -> Double {
  total_pressure_plane(rho, G, h_c, area)
}

///|
/// 计算压力中心深度 (矩形平面，垂直放置)
/// y_p = y_c + I_c / (y_c * A)
/// 参数: y_c - 形心深度 (m), i_c - 对形心轴的惯性矩 (m⁴), area - 面积 (m²)
/// 返回: 压力中心深度 (m)
pub fn pressure_center_depth(
  y_c : Double,
  i_c : Double,
  area : Double,
) -> Double {
  if area <= 0.0 || y_c <= 0.0 {
    0.0
  } else {
    y_c + i_c / (y_c * area)
  }
}

///|
/// 矩形平面压力中心深度
/// 参数: h_c - 形心深度 (m), b - 宽度 (m), h - 高度 (m)
pub fn pressure_center_rect(h_c : Double, b : Double, h : Double) -> Double {
  let i_c = b * cube(h) / 12.0
  let area = b * h
  pressure_center_depth(h_c, i_c, area)
}

// ==================== 连续性方程 ====================

///|
/// 连续性方程 - 体积流量守恒
/// Q = A₁v₁ = A₂v₂
/// 参数: area1 - 截面1面积 (m²), v1 - 截面1流速 (m/s)
/// 返回: 体积流量 (m³/s)
pub fn volume_flow_rate(area : Double, velocity : Double) -> Double {
  area * velocity
}

///|
/// 连续性方程 - 计算流速
/// v = Q / A
pub fn velocity_from_flow_rate(flow_rate : Double, area : Double) -> Double {
  if area <= 0.0 {
    0.0
  } else {
    flow_rate / area
  }
}

///|
/// 连续性方程 - 计算截面面积
/// A = Q / v
pub fn area_from_flow_rate(flow_rate : Double, velocity : Double) -> Double {
  if velocity <= 0.0 {
    0.0
  } else {
    flow_rate / velocity
  }
}

///|
/// 质量流量
/// ṁ = ρQ = ρAv
pub fn mass_flow_rate(rho : Double, flow_rate : Double) -> Double {
  rho * flow_rate
}

///|
/// 质量流量 (从面积和流速)
pub fn mass_flow_rate_av(
  rho : Double,
  area : Double,
  velocity : Double,
) -> Double {
  rho * area * velocity
}

// ==================== 伯努利方程 ====================

///|
/// 伯努利方程 - 总压头
/// H = z + p/(ρg) + v²/(2g)
/// 参数: z - 位置水头 (m), p - 压强 (Pa), rho - 密度 (kg/m³), g - 重力加速度 (m/s²), v - 流速 (m/s)
/// 返回: 总压头 (m)
pub fn bernoulli_total_head(
  z : Double,
  p : Double,
  rho : Double,
  g : Double,
  v : Double,
) -> Double {
  z + p / (rho * g) + square(v) / (2.0 * g)
}

///|
/// 伯努利方程 - 总压头 (使用标准重力)
pub fn bernoulli_total_head_std(
  z : Double,
  p : Double,
  rho : Double,
  v : Double,
) -> Double {
  bernoulli_total_head(z, p, rho, G, v)
}

///|
/// 伯努利方程 - 压强水头
/// h_p = p / (ρg)
pub fn pressure_head(p : Double, rho : Double, g : Double) -> Double {
  if rho <= 0.0 || g <= 0.0 {
    0.0
  } else {
    p / (rho * g)
  }
}

///|
/// 伯努利方程 - 速度水头
/// h_v = v² / (2g)
pub fn velocity_head(v : Double, g : Double) -> Double {
  if g <= 0.0 {
    0.0
  } else {
    square(v) / (2.0 * g)
  }
}

///|
/// 伯努利方程 - 两点间压差
/// p₁ + ½ρv₁² + ρgz₁ = p₂ + ½ρv₂² + ρgz₂
/// 返回: p₁ - p₂
pub fn bernoulli_pressure_drop(
  v1 : Double,
  v2 : Double,
  z1 : Double,
  z2 : Double,
  rho : Double,
  g : Double,
) -> Double {
  rho * g * (z2 - z1) + 0.5 * rho * (square(v2) - square(v1))
}

///|
/// 文丘里流量计 - 计算流量
/// Q = C * A₂ * sqrt(2gΔh / (1 - (A₂/A₁)²))
/// 参数: c - 流量系数, a1 - 上游面积, a2 - 喉部面积, dh - 压差水头 (m), g - 重力加速度
pub fn venturi_flow_rate(
  c : Double,
  a1 : Double,
  a2 : Double,
  dh : Double,
  g : Double,
) -> Double {
  if a1 <= 0.0 || a2 <= 0.0 || dh <= 0.0 {
    0.0
  } else {
    let beta = a2 / a1
    let denom = 1.0 - square(beta)
    if denom <= 0.0 {
      0.0
    } else {
      c * a2 * (2.0 * g * dh / denom).sqrt()
    }
  }
}

///|
/// 孔板流量计 - 计算流量
/// Q = C * A₀ * sqrt(2gΔh)
/// 参数: c - 流量系数, a0 - 孔口面积, dh - 压差水头 (m), g - 重力加速度
pub fn orifice_flow_rate(
  c : Double,
  a0 : Double,
  dh : Double,
  g : Double,
) -> Double {
  if a0 <= 0.0 || dh <= 0.0 {
    0.0
  } else {
    c * a0 * (2.0 * g * dh).sqrt()
  }
}

// ==================== 动量方程 ====================

///|
/// 动量方程 - 作用力
/// F = ṁ(v₂ - v₁) = ρQ(v₂ - v₁)
/// 参数: rho - 密度 (kg/m³), q - 流量 (m³/s), v1 - 入口流速 (m/s), v2 - 出口流速 (m/s)
/// 返回: 作用力 (N)
pub fn momentum_force(
  rho : Double,
  q : Double,
  v1 : Double,
  v2 : Double,
) -> Double {
  rho * q * (v2 - v1)
}

///|
/// 射流冲击力 (垂直冲击平板)
/// F = ρQv = ρAv²
pub fn jet_force(rho : Double, area : Double, velocity : Double) -> Double {
  rho * area * square(velocity)
}

///|
/// 射流冲击力 (倾斜冲击平板，角度为法向)
/// F = ρQv * sin(θ)
pub fn jet_force_angled(
  rho : Double,
  area : Double,
  velocity : Double,
  angle_rad : Double,
) -> Double {
  // 使用近似 sin
  let sin_theta = sin_approx(angle_rad)
  rho * area * square(velocity) * sin_theta
}

///|
/// 正弦函数近似
fn sin_approx(x : Double) -> Double {
  let mut a = x
  while a > PI {
    a = a - 2.0 * PI
  }
  while a < -PI {
    a = a + 2.0 * PI
  }
  let x2 = a * a
  let x3 = x2 * a
  let x5 = x3 * x2
  let x7 = x5 * x2
  a - x3 / 6.0 + x5 / 120.0 - x7 / 5040.0
}

// ==================== 雷诺数 ====================

///|
/// 雷诺数 (圆管)
/// Re = ρvD / μ = vD / ν
/// 参数: rho - 密度 (kg/m³), v - 流速 (m/s), d - 管径 (m), mu - 动力粘度 (Pa·s)
/// 返回: 雷诺数 (无量纲)
pub fn reynolds_number(
  rho : Double,
  v : Double,
  d : Double,
  mu : Double,
) -> Double {
  if mu <= 0.0 {
    0.0
  } else {
    rho * v * d / mu
  }
}

///|
/// 雷诺数 (使用运动粘度)
/// Re = vD / ν
pub fn reynolds_number_kinematic(v : Double, d : Double, nu : Double) -> Double {
  if nu <= 0.0 {
    0.0
  } else {
    v * d / nu
  }
}

///|
/// 雷诺数 (非圆管，使用水力直径)
/// Re = ρvD_h / μ
pub fn reynolds_number_hydraulic_diameter(
  rho : Double,
  v : Double,
  d_h : Double,
  mu : Double,
) -> Double {
  reynolds_number(rho, v, d_h, mu)
}

///|
/// 水力直径
/// D_h = 4A / P
/// 参数: area - 过流面积 (m²), perimeter - 湿周 (m)
pub fn hydraulic_diameter(area : Double, perimeter : Double) -> Double {
  if perimeter <= 0.0 {
    0.0
  } else {
    4.0 * area / perimeter
  }
}

///|
/// 判断流动状态 (层流/湍流)
/// 返回: true 为湍流，false 为层流
pub fn is_turbulent(re : Double, re_critical : Double) -> Bool {
  re > re_critical
}

///|
/// 判断圆管流动状态
pub fn is_turbulent_pipe(re : Double) -> Bool {
  is_turbulent(re, RE_CRITICAL_PIPE)
}

// ==================== 管道流动 - 层流 ====================

///|
/// 圆管层流速度分布 (抛物线)
/// u(r) = (Δp/(4μL)) * (R² - r²)
/// 参数: dp - 压降 (Pa), mu - 动力粘度 (Pa·s), l - 管长 (m), r - 半径 (m), r_pipe - 管半径 (m)
/// 返回: 流速 (m/s)
pub fn laminar_velocity_profile(
  dp : Double,
  mu : Double,
  l : Double,
  r : Double,
  r_pipe : Double,
) -> Double {
  if mu <= 0.0 || l <= 0.0 || r > r_pipe {
    0.0
  } else {
    dp / (4.0 * mu * l) * (square(r_pipe) - square(r))
  }
}

///|
/// 圆管层流最大流速 (管中心)
pub fn laminar_max_velocity(
  dp : Double,
  mu : Double,
  l : Double,
  r : Double,
) -> Double {
  if mu <= 0.0 || l <= 0.0 {
    0.0
  } else {
    dp * square(r) / (4.0 * mu * l)
  }
}

///|
/// 圆管层流平均流速
/// v_avg = v_max / 2 = ΔpR² / (8μL)
pub fn laminar_mean_velocity(
  dp : Double,
  mu : Double,
  l : Double,
  r : Double,
) -> Double {
  laminar_max_velocity(dp, mu, l, r) / 2.0
}

///|
/// 圆管层流流量 (Hagen-Poiseuille)
/// Q = πΔpR⁴ / (8μL)
pub fn laminar_flow_rate(
  dp : Double,
  mu : Double,
  l : Double,
  r : Double,
) -> Double {
  if mu <= 0.0 || l <= 0.0 {
    0.0
  } else {
    PI * dp * pow4(r) / (8.0 * mu * l)
  }
}

///|
/// 圆管层流压降
/// Δp = 32μLv / D²
pub fn laminar_pressure_drop(
  mu : Double,
  l : Double,
  v : Double,
  d : Double,
) -> Double {
  if d <= 0.0 {
    0.0
  } else {
    32.0 * mu * l * v / square(d)
  }
}

///|
/// 圆管层流摩擦系数
/// f = 64 / Re
pub fn laminar_friction_factor(re : Double) -> Double {
  if abs(re) < EPSILON {
    0.0
  } else {
    64.0 / re
  }
}

// ==================== 管道流动 - 湍流 ====================

///|
/// 达西-魏斯巴赫方程 - 摩擦损失水头
/// h_f = f * (L/D) * (v²/(2g))
/// 参数: f - 摩擦系数, l - 管长 (m), d - 管径 (m), v - 流速 (m/s), g - 重力加速度
pub fn darcy_weisbach_head_loss(
  f : Double,
  l : Double,
  d : Double,
  v : Double,
  g : Double,
) -> Double {
  if d <= 0.0 || g <= 0.0 {
    0.0
  } else {
    f * (l / d) * square(v) / (2.0 * g)
  }
}

///|
/// 达西-魏斯巴赫方程 - 摩擦损失水头 (使用标准重力)
pub fn darcy_weisbach_head_loss_std(
  f : Double,
  l : Double,
  d : Double,
  v : Double,
) -> Double {
  darcy_weisbach_head_loss(f, l, d, v, G)
}

///|
/// 达西-魏斯巴赫方程 - 压降
/// Δp = f * (L/D) * (ρv²/2)
pub fn darcy_weisbach_pressure_drop(
  f : Double,
  l : Double,
  d : Double,
  rho : Double,
  v : Double,
) -> Double {
  if d <= 0.0 {
    0.0
  } else {
    f * (l / d) * rho * square(v) / 2.0
  }
}

///|
/// 科尔布鲁克方程 - 摩擦系数 (隐式，需迭代求解)
/// 1/√f = -2 log₁₀(ε/(3.7D) + 2.51/(Re√f))
/// 这里提供初始估计值
pub fn colebrook_friction_factor_estimate(
  re : Double,
  epsilon : Double,
  d : Double,
) -> Double {
  if re <= 0.0 || d <= 0.0 {
    return 0.0
  }
  let relative_roughness = epsilon / d
  // Swamee-Jain 近似公式
  let term1 = relative_roughness / 3.7
  let term2 = 5.74 / pow_int(re, 9)
  let f = 0.25 / square(log10_approx(term1 + term2))
  f
}

///|
/// 莫迪图 - 光滑管摩擦系数 (Blasius公式，Re < 10⁵)
/// f = 0.316 / Re^0.25
pub fn blasius_friction_factor(re : Double) -> Double {
  if re <= 0.0 {
    0.0
  } else {
    0.316 / pow_int(re, 4).sqrt().sqrt() // Re^0.25
  }
}

///|
/// 完全粗糙管摩擦系数 (Re 很大时)
/// f = 0.25 / (log₁₀(ε/(3.7D)))²
pub fn fully_rough_friction_factor(epsilon : Double, d : Double) -> Double {
  if d <= 0.0 || epsilon <= 0.0 {
    0.0
  } else {
    let term = log10_approx(epsilon / (3.7 * d))
    0.25 / square(term)
  }
}

// ==================== 局部损失 ====================

///|
/// 局部损失水头
/// h_l = K * v² / (2g)
/// 参数: k - 局部阻力系数, v - 流速 (m/s), g - 重力加速度
pub fn local_head_loss(k : Double, v : Double, g : Double) -> Double {
  if g <= 0.0 {
    0.0
  } else {
    k * square(v) / (2.0 * g)
  }
}

///|
/// 局部损失水头 (使用标准重力)
pub fn local_head_loss_std(k : Double, v : Double) -> Double {
  local_head_loss(k, v, G)
}

///|
/// 突然扩大损失系数
/// K = (1 - A₁/A₂)²
pub fn sudden_expansion_loss_coefficient(a1 : Double, a2 : Double) -> Double {
  if a2 <= 0.0 {
    0.0
  } else {
    square(1.0 - a1 / a2)
  }
}

///|
/// 突然缩小损失系数 (经验公式)
/// K ≈ 0.5 * (1 - A₂/A₁)
pub fn sudden_contraction_loss_coefficient(a1 : Double, a2 : Double) -> Double {
  if a1 <= 0.0 {
    0.0
  } else {
    0.5 * (1.0 - a2 / a1)
  }
}

///|
/// 总水头损失 (摩擦 + 局部)
pub fn total_head_loss(h_friction : Double, h_local : Double) -> Double {
  h_friction + h_local
}

// ==================== 明渠流 ====================

///|
/// 曼宁公式 - 流量
/// Q = (1/n) * A * R^(2/3) * S^(1/2)
/// 参数: n - 曼宁粗糙度系数, area - 过流面积 (m²), r - 水力半径 (m), s - 底坡
pub fn manning_flow_rate(
  n : Double,
  area : Double,
  r : Double,
  s : Double,
) -> Double {
  if n <= 0.0 || area <= 0.0 || r <= 0.0 || s < 0.0 {
    0.0
  } else {
    // R^(2/3) = (R^2)^(1/3) = (R^2).cbrt()
    let r_23 = cbrt(pow_int(r, 2))
    1.0 / n * area * r_23 * s.sqrt()
  }
}

///|
/// 曼宁公式 - 流速
/// v = (1/n) * R^(2/3) * S^(1/2)
pub fn manning_velocity(n : Double, r : Double, s : Double) -> Double {
  if n <= 0.0 || r <= 0.0 || s < 0.0 {
    0.0
  } else {
    let r_23 = cbrt(pow_int(r, 2))
    1.0 / n * r_23 * s.sqrt()
  }
}

///|
/// 水力半径
/// R = A / P
/// 参数: area - 过流面积 (m²), perimeter - 湿周 (m)
pub fn hydraulic_radius(area : Double, perimeter : Double) -> Double {
  if perimeter <= 0.0 {
    0.0
  } else {
    area / perimeter
  }
}

///|
/// 矩形明渠水力半径
pub fn hydraulic_radius_rect(b : Double, h : Double) -> Double {
  if b <= 0.0 || h <= 0.0 {
    0.0
  } else {
    let area = b * h
    let perimeter = b + 2.0 * h
    hydraulic_radius(area, perimeter)
  }
}

///|
/// 梯形明渠水力半径
/// 参数: b - 底宽 (m), h - 水深 (m), m - 边坡系数 (水平/垂直)
pub fn hydraulic_radius_trapezoid(b : Double, h : Double, m : Double) -> Double {
  if b <= 0.0 || h <= 0.0 {
    0.0
  } else {
    let area = (b + m * h) * h
    let perimeter = b + 2.0 * h * (1.0 + square(m)).sqrt()
    hydraulic_radius(area, perimeter)
  }
}

///|
/// 临界水深 (矩形明渠)
/// h_c = (q²/g)^(1/3)
/// 参数: q - 单宽流量 (m²/s), g - 重力加速度
pub fn critical_depth_rect(q : Double, g : Double) -> Double {
  if g <= 0.0 || q <= 0.0 {
    0.0
  } else {
    cbrt(square(q) / g) // (q²/g)^(1/3)
  }
}

///|
/// 弗劳德数
/// Fr = v / sqrt(gh)
/// 参数: v - 流速 (m/s), g - 重力加速度, h - 特征水深 (m)
pub fn froude_number(v : Double, g : Double, h : Double) -> Double {
  if g <= 0.0 || h <= 0.0 {
    0.0
  } else {
    v / (g * h).sqrt()
  }
}

///|
/// 判断流动状态 (Fr < 1 缓流, Fr = 1 临界, Fr > 1 急流)
/// 返回: "subcritical", "critical", "supercritical"
pub fn flow_regime(fr : Double) -> String {
  if abs(fr - 1.0) < 0.01 {
    "critical"
  } else if fr < 1.0 {
    "subcritical"
  } else {
    "supercritical"
  }
}

// ==================== 堰流 ====================

///|
/// 矩形薄壁堰流量 (Francis公式)
/// Q = C * b * H^(3/2)
/// 参数: c - 流量系数 (约 1.84), b - 堰宽 (m), h - 堰上水头 (m)
pub fn weir_rectangular_flow(c : Double, b : Double, h : Double) -> Double {
  if h <= 0.0 || b <= 0.0 {
    0.0
  } else {
    c * b * h.sqrt() * h // H^(3/2) = H * sqrt(H)
  }
}

///|
/// 三角形薄壁堰流量
/// Q = C * tan(θ/2) * H^(5/2)
/// 参数: c - 流量系数 (约 1.4), theta_deg - 堰口角度 (度), h - 堰上水头 (m)
pub fn weir_triangular_flow(
  c : Double,
  theta_deg : Double,
  h : Double,
) -> Double {
  if h <= 0.0 {
    0.0
  } else {
    // tan(θ/2) 近似
    let theta_rad = theta_deg * PI / 180.0
    let tan_half = sin_approx(theta_rad / 2.0) / cos_approx(theta_rad / 2.0)
    c * tan_half * square(h) * h.sqrt() // H^(5/2) = H^2 * sqrt(H)
  }
}

///|
/// 余弦函数近似
fn cos_approx(x : Double) -> Double {
  let mut a = x
  while a > PI {
    a = a - 2.0 * PI
  }
  while a < -PI {
    a = a + 2.0 * PI
  }
  let x2 = a * a
  let x4 = x2 * x2
  let x6 = x4 * x2
  1.0 - x2 / 2.0 + x4 / 24.0 - x6 / 720.0
}

// ==================== 绕流问题 ====================

///|
/// 阻力
/// F_D = ½ * ρ * v² * A * C_D
/// 参数: rho - 密度 (kg/m³), v - 流速 (m/s), area - 迎风面积 (m²), cd - 阻力系数
/// 返回: 阻力 (N)
pub fn drag_force(
  rho : Double,
  v : Double,
  area : Double,
  cd : Double,
) -> Double {
  0.5 * rho * square(v) * area * cd
}

///|
/// 升力
/// F_L = ½ * ρ * v² * A * C_L
pub fn lift_force(
  rho : Double,
  v : Double,
  area : Double,
  cl : Double,
) -> Double {
  0.5 * rho * square(v) * area * cl
}

///|
/// 球体阻力系数 (Stokes定律，Re < 1)
/// C_D = 24 / Re
pub fn sphere_drag_coefficient_stokes(re : Double) -> Double {
  if abs(re) < EPSILON {
    0.0
  } else {
    24.0 / re
  }
}

///|
/// 球体阻力 (Stokes定律)
/// F_D = 3πμDv
pub fn sphere_drag_stokes(mu : Double, d : Double, v : Double) -> Double {
  3.0 * PI * mu * d * v
}

///|
/// 圆柱绕流阻力系数 (经验值，Re ~ 10³-10⁵)
pub fn cylinder_drag_coefficient(re : Double) -> Double {
  if re < 1.0 {
    10.0 / re.sqrt()
  } else if re < 1000.0 {
    1.0 + 10.0 / re.sqrt()
  } else if re < 100000.0 {
    1.2
  } else {
    0.3
  }
}

// ==================== 量纲分析 ====================

///|
/// 计算雷诺数 (量纲分析)
pub fn reynolds_dimensionless(
  rho : Double,
  v : Double,
  l : Double,
  mu : Double,
) -> Double {
  reynolds_number(rho, v, l, mu)
}

///|
/// 计算弗劳德数 (量纲分析)
pub fn froude_dimensionless(v : Double, g : Double, l : Double) -> Double {
  froude_number(v, g, l)
}

///|
/// 计算韦伯数
/// We = ρv²L / σ
/// 参数: rho - 密度, v - 流速, l - 特征长度, sigma - 表面张力
pub fn weber_number(
  rho : Double,
  v : Double,
  l : Double,
  sigma : Double,
) -> Double {
  if sigma <= 0.0 {
    0.0
  } else {
    rho * square(v) * l / sigma
  }
}

///|
/// 计算马赫数
/// Ma = v / c
/// 参数: v - 流速, c - 声速
pub fn mach_number(v : Double, c : Double) -> Double {
  if c <= 0.0 {
    0.0
  } else {
    v / c
  }
}

///|
/// 计算欧拉数
/// Eu = Δp / (ρv²)
pub fn euler_number(dp : Double, rho : Double, v : Double) -> Double {
  if rho <= 0.0 || abs(v) < EPSILON {
    0.0
  } else {
    dp / (rho * square(v))
  }
}

// ==================== 相似理论 ====================

///|
/// 检查雷诺数相似
pub fn reynolds_similarity(re_model : Double, re_prototype : Double) -> Bool {
  approx_equal(re_model, re_prototype, 0.01)
}

///|
/// 检查弗劳德数相似
pub fn froude_similarity(fr_model : Double, fr_prototype : Double) -> Bool {
  approx_equal(fr_model, fr_prototype, 0.01)
}

// ==================== 边界层 ====================

///|
/// 平板边界层厚度 (层流，Blasius解)
/// δ = 5.0 * x / sqrt(Re_x)
/// 参数: x - 距前缘距离 (m), re_x - 当地雷诺数
pub fn boundary_layer_thickness_laminar(x : Double, re_x : Double) -> Double {
  if re_x <= 0.0 {
    0.0
  } else {
    5.0 * x / re_x.sqrt()
  }
}

///|
/// 平板边界层厚度 (湍流，1/7次方律)
/// δ = 0.37 * x / Re_x^(1/5)
pub fn boundary_layer_thickness_turbulent(x : Double, re_x : Double) -> Double {
  if re_x <= 0.0 {
    0.0
  } else {
    // Re_x^(1/5) 近似计算
    let re_15 = cbrt(cbrt(cbrt(pow_int(re_x, 1)))) // 近似 Re^(1/5)
    0.37 * x / re_15
  }
}

///|
/// 平板摩擦阻力系数 (层流)
/// C_f = 0.664 / sqrt(Re_x)
pub fn flat_plate_friction_coefficient_laminar(re_x : Double) -> Double {
  if re_x <= 0.0 {
    0.0
  } else {
    0.664 / re_x.sqrt()
  }
}

///|
/// 平板摩擦阻力系数 (湍流)
/// C_f = 0.0592 / Re_x^(1/5)
pub fn flat_plate_friction_coefficient_turbulent(re_x : Double) -> Double {
  if re_x <= 0.0 {
    0.0
  } else {
    let re_15 = cbrt(cbrt(cbrt(pow_int(re_x, 1)))) // 近似 Re^(1/5)
    0.0592 / re_15
  }
}
