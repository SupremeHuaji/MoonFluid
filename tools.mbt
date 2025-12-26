// MoonFluid - 工具函数

// ==================== 基本数学运算 ====================

///|
/// 计算平方
pub fn square(x : Double) -> Double {
  x * x
}

///|
/// 计算立方
pub fn cube(x : Double) -> Double {
  x * x * x
}

///|
/// 计算四次方
pub fn pow4(x : Double) -> Double {
  let x2 = x * x
  x2 * x2
}

///|
/// 计算绝对值
pub fn abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
/// 计算最大值
pub fn max(a : Double, b : Double) -> Double {
  if a > b {
    a
  } else {
    b
  }
}

///|
/// 计算最小值
pub fn min(a : Double, b : Double) -> Double {
  if a < b {
    a
  } else {
    b
  }
}

///|
/// 安全除法，分母为零时返回默认值
pub fn safe_div(num : Double, den : Double, fallback : Double) -> Double {
  if abs(den) <= EPSILON {
    fallback
  } else {
    num / den
  }
}

///|
/// 判断两个浮点数是否近似相等
pub fn approx_equal(a : Double, b : Double, tolerance : Double) -> Bool {
  abs(a - b) < tolerance
}

///|
/// 判断浮点数是否接近零
pub fn is_zero(x : Double) -> Bool {
  abs(x) < EPSILON
}

///|
/// 符号函数
pub fn sign(x : Double) -> Double {
  if x > EPSILON {
    1.0
  } else if x < -EPSILON {
    -1.0
  } else {
    0.0
  }
}

///|
/// 限制值在范围内
pub fn clamp(value : Double, min_val : Double, max_val : Double) -> Double {
  if value < min_val {
    min_val
  } else if value > max_val {
    max_val
  } else {
    value
  }
}

///|
/// 线性插值
pub fn lerp(a : Double, b : Double, t : Double) -> Double {
  a + (b - a) * t
}

// ==================== 幂函数 ====================

///|
/// 计算 x 的 n 次方 (n 为整数)
pub fn pow_int(x : Double, n : Int) -> Double {
  if n == 0 {
    1.0
  } else if n < 0 {
    1.0 / pow_int(x, -n)
  } else {
    let mut result = 1.0
    let mut base = x
    let mut exp = n
    while exp > 0 {
      if exp % 2 == 1 {
        result = result * base
      }
      base = base * base
      exp = exp / 2
    }
    result
  }
}

///|
/// 计算 x 的 1/n 次方 (n 次根)
pub fn root(x : Double, n : Int) -> Double {
  if n <= 0 {
    0.0
  } else if n == 1 {
    x
  } else if n == 2 {
    x.sqrt()
    // 使用牛顿法: x^(1/n) = y, 其中 y^n = x
    // 迭代: y_new = y - (y^n - x) / (n * y^(n-1))
  } else if x <= 0.0 {
    0.0
  } else {
    let mut y = x
    let mut i = 0
    while i < 20 {
      let yn = pow_int(y, n)
      let dydn = n.to_double() * pow_int(y, n - 1)
      let y_new = y - (yn - x) / dydn
      if abs(y_new - y) < EPSILON {
        return y_new
      }
      y = y_new
      i = i + 1
    }
    y
  }
}

// ==================== 几何计算 ====================

///|
/// 圆形面积
pub fn circle_area(r : Double) -> Double {
  PI * square(r)
}

///|
/// 圆形周长
pub fn circle_perimeter(r : Double) -> Double {
  2.0 * PI * r
}

///|
/// 矩形面积
pub fn rect_area(width : Double, height : Double) -> Double {
  width * height
}

///|
/// 梯形面积
pub fn trapezoid_area(a : Double, b : Double, h : Double) -> Double {
  (a + b) * h / 2.0
}

///|
/// 圆环面积
pub fn ring_area(r_outer : Double, r_inner : Double) -> Double {
  PI * (square(r_outer) - square(r_inner))
}

// ==================== 统计工具 ====================

///|
/// 数组求和
pub fn array_sum(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  for i = 0; i < arr.length(); i = i + 1 {
    sum = sum + arr[i]
  }
  sum
}

///|
/// 数组平均值
pub fn array_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    0.0
  } else {
    array_sum(arr) / arr.length().to_double()
  }
}

///|
/// 数组最大值
pub fn array_max(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut max_val = arr[0]
  for i = 1; i < arr.length(); i = i + 1 {
    if arr[i] > max_val {
      max_val = arr[i]
    }
  }
  max_val
}

///|
/// 数组最小值
pub fn array_min(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut min_val = arr[0]
  for i = 1; i < arr.length(); i = i + 1 {
    if arr[i] < min_val {
      min_val = arr[i]
    }
  }
  min_val
}

// ==================== 数值方法 ====================

///|
/// 二分法求根 (f(x) = 0)
/// 函数 f 在 [a, b] 上连续，且 f(a)*f(b) < 0
pub fn bisection(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  tolerance : Double,
) -> Double {
  let mut left = a
  let mut right = b
  let mut mid = (left + right) / 2.0
  let mut i = 0
  while abs(right - left) > tolerance && i < 100 {
    mid = (left + right) / 2.0
    let f_mid = f(mid)
    let f_left = f(left)
    if f_mid * f_left < 0.0 {
      right = mid
    } else {
      left = mid
    }
    i = i + 1
  }
  mid
}

///|
/// 牛顿法求根
/// f: 函数, df: 导数, x0: 初始值
pub fn newton_method(
  f : (Double) -> Double,
  df : (Double) -> Double,
  x0 : Double,
  tolerance : Double,
) -> Double {
  let mut x = x0
  let mut i = 0
  while i < 50 {
    let fx = f(x)
    if abs(fx) < tolerance {
      return x
    }
    let dfx = df(x)
    if abs(dfx) < EPSILON {
      return x
    }
    x = x - fx / dfx
    i = i + 1
  }
  x
}

// ==================== 插值 ====================

///|
/// 自然对数近似 (使用级数展开)
pub fn ln_approx(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x == 1.0 {
    return 0.0
  }
  // 对于 x > 1，使用 ln(x) = -ln(1/x)
  if x > 1.0 {
    return -ln_approx(1.0 / x)
  }
  // 对于 0 < x < 1，使用级数: ln(1+u) = u - u²/2 + u³/3 - ...
  let u = x - 1.0
  let mut sum = u
  let mut term = u
  let mut i = 2
  while i < 50 {
    term = term * -u
    sum = sum + term / i.to_double()
    if abs(term) < EPSILON {
      break
    }
    i = i + 1
  }
  sum
}

///|
/// 以10为底的对数
pub fn log10_approx(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  ln_approx(x) / 2.302585092994046 // ln(10)
}

///|
/// 立方根
pub fn cbrt(x : Double) -> Double {
  if x == 0.0 {
    0.0
  } else if x < 0.0 {
    -cbrt(-x)
  } else {
    // 使用牛顿法: y = x^(1/3), y^3 = x
    // 迭代: y_new = y - (y^3 - x) / (3*y^2)
    let mut y = x
    let mut i = 0
    while i < 20 {
      let y3 = y * y * y
      let y2 = y * y
      let y_new = y - (y3 - x) / (3.0 * y2)
      if abs(y_new - y) < EPSILON {
        return y_new
      }
      y = y_new
      i = i + 1
    }
    y
  }
}

///|
/// 线性插值表查找
/// x_values: 自变量数组 (需单调递增)
/// y_values: 因变量数组
/// x: 待插值的自变量
pub fn linear_interpolate(
  x_values : Array[Double],
  y_values : Array[Double],
  x : Double,
) -> Double {
  let n = x_values.length()
  if n == 0 {
    return 0.0
  }
  if n == 1 {
    return y_values[0]
  }
  if x <= x_values[0] {
    return y_values[0]
  }
  if x >= x_values[n - 1] {
    return y_values[n - 1]
  }
  // 找到 x 所在区间
  let mut i = 0
  while i < n - 1 {
    if x >= x_values[i] && x <= x_values[i + 1] {
      let t = (x - x_values[i]) / (x_values[i + 1] - x_values[i])
      return lerp(y_values[i], y_values[i + 1], t)
    }
    i = i + 1
  }
  y_values[n - 1]
}
