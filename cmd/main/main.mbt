///|
fn main {
  // 示例：计算管道流量
  // 圆管层流: Δp=5000 Pa, μ=0.001 Pa·s, L=20 m, R=0.05 m
  let q = @lib.laminar_flow_rate(5000.0, 0.001, 20.0, 0.05)
  println("层流流量: \{q} m³/s")

  // 示例：计算雷诺数
  let re = @lib.reynolds_number(1000.0, 2.0, 0.1, 0.001)
  println("雷诺数: \{re}")

  // 示例：曼宁公式计算明渠流量
  let q_manning = @lib.manning_flow_rate(0.013, 5.0, 1.0, 0.001)
  println("明渠流量: \{q_manning} m³/s")
}
