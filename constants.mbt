// MoonFluid - 流体力学常量定义

// ==================== 数学常量 ====================

///|
/// 圆周率
pub const PI : Double = 3.141592653589793

///|
/// 自然对数底数
pub const E_CONST : Double = 2.718281828459045

///|
/// 浮点数比较精度
pub const EPSILON : Double = 1.0e-10

// ==================== 物理常量 ====================

///|
/// 标准重力加速度 (m/s²)
pub const G : Double = 9.80665

///|
/// 标准大气压 (Pa)
pub const ATMOSPHERIC_PRESSURE : Double = 101325.0

///|
/// 水的密度 (kg/m³) @ 4°C
pub const RHO_WATER : Double = 1000.0

///|
/// 空气密度 (kg/m³) @ 标准状态 (20°C, 1 atm)
pub const RHO_AIR : Double = 1.204

///|
/// 海水密度 (kg/m³) @ 15°C
pub const RHO_SEAWATER : Double = 1025.0

///|
/// 水银密度 (kg/m³)
pub const RHO_MERCURY : Double = 13595.0

///|
/// 汽油密度 (kg/m³)
pub const RHO_GASOLINE : Double = 750.0

///|
/// 柴油密度 (kg/m³)
pub const RHO_DIESEL : Double = 850.0

// ==================== 流体物性参数 ====================

///|
/// 水的动力粘度 (Pa·s) @ 20°C
pub const MU_WATER_20C : Double = 0.001002

///|
/// 水的运动粘度 (m²/s) @ 20°C
pub const NU_WATER_20C : Double = 1.004e-6

///|
/// 空气的动力粘度 (Pa·s) @ 20°C
pub const MU_AIR_20C : Double = 1.81e-5

///|
/// 空气的运动粘度 (m²/s) @ 20°C
pub const NU_AIR_20C : Double = 1.51e-5

///|
/// 水的表面张力 (N/m) @ 20°C
pub const SIGMA_WATER : Double = 0.0728

///|
/// 水的体积模量 (Pa) @ 20°C
pub const K_WATER : Double = 2.2e9

// ==================== 雷诺数临界值 ====================

///|
/// 圆管层流临界雷诺数
pub const RE_CRITICAL_PIPE : Double = 2300.0

///|
/// 平板边界层层流临界雷诺数
pub const RE_CRITICAL_PLATE : Double = 5.0e5

///|
/// 明渠流临界雷诺数
pub const RE_CRITICAL_CHANNEL : Double = 500.0

// ==================== 粗糙度系数 ====================

///|
/// 曼宁粗糙度系数 - 光滑混凝土
pub const N_MANNING_SMOOTH_CONCRETE : Double = 0.012

///|
/// 曼宁粗糙度系数 - 普通混凝土
pub const N_MANNING_CONCRETE : Double = 0.013

///|
/// 曼宁粗糙度系数 - 砖砌
pub const N_MANNING_BRICK : Double = 0.015

///|
/// 曼宁粗糙度系数 - 土渠
pub const N_MANNING_EARTH : Double = 0.025

///|
/// 曼宁粗糙度系数 - 天然河道
pub const N_MANNING_NATURAL : Double = 0.035

///|
/// 曼宁粗糙度系数 - 杂草河道
pub const N_MANNING_WEEDY : Double = 0.040

// ==================== 管道粗糙度 (绝对粗糙度, m) ====================

///|
/// 光滑管 (玻璃、塑料)
pub const EPSILON_SMOOTH : Double = 0.0

///|
/// 新钢管
pub const EPSILON_STEEL_NEW : Double = 0.000045

///|
/// 旧钢管
pub const EPSILON_STEEL_OLD : Double = 0.0002

///|
/// 铸铁管
pub const EPSILON_CAST_IRON : Double = 0.00026

///|
/// 混凝土管
pub const EPSILON_CONCRETE : Double = 0.0003

///|
/// 旧铸铁管
pub const EPSILON_CAST_IRON_OLD : Double = 0.0015

// ==================== 局部阻力系数 ====================

///|
/// 90° 弯头 (标准)
pub const K_BEND_90 : Double = 0.9

///|
/// 45° 弯头
pub const K_BEND_45 : Double = 0.4

///|
/// 突然扩大 (A2/A1 → 0)
pub const K_SUDDEN_EXPANSION_MAX : Double = 1.0

///|
/// 突然缩小 (A2/A1 → 0)
pub const K_SUDDEN_CONTRACTION_MAX : Double = 0.5

///|
/// 阀门 (全开)
pub const K_VALVE_OPEN : Double = 0.2

///|
/// 阀门 (半开)
pub const K_VALVE_HALF : Double = 4.0

///|
/// 入口 (锐缘)
pub const K_INLET_SHARP : Double = 0.5

///|
/// 入口 (圆滑)
pub const K_INLET_ROUNDED : Double = 0.05

///|
/// 出口
pub const K_OUTLET : Double = 1.0

// ==================== 单位换算 ====================

///|
/// 1 bar = 10⁵ Pa
pub const BAR_TO_PA : Double = 100000.0

///|
/// 1 atm = 101325 Pa
pub const ATM_TO_PA : Double = 101325.0

///|
/// 1 psi = 6894.76 Pa
pub const PSI_TO_PA : Double = 6894.76

///|
/// 1 m³/s = 1000 L/s
pub const M3S_TO_LS : Double = 1000.0

///|
/// 1 L/s = 0.001 m³/s
pub const LS_TO_M3S : Double = 0.001

///|
/// 1 m³/h = 1/3600 m³/s
pub const M3H_TO_M3S : Double = 1.0 / 3600.0

// ==================== 经验系数 ====================

///|
/// 文丘里流量计系数 (典型值)
pub const C_VENTURI : Double = 0.98

///|
/// 孔板流量计系数 (典型值)
pub const C_ORIFICE : Double = 0.61

///|
/// 皮托管系数 (典型值)
pub const C_PITOT : Double = 1.0
