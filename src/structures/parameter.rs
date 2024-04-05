use std::ops;

#[derive(Clone, Debug, PartialEq)]
pub enum Parameter {
    Abstract(String),
    Concrete(f64),
}

impl Parameter {
    pub fn zero() -> Self {
        return Self::Concrete(0.0);
    }
    pub fn is_zero(&self) -> bool {
        match self {
            Self::Abstract(_) => false,
            Self::Concrete(x) => x.abs() < 1e-6,
        }
    }
    pub fn flip_sign(&mut self) {
        match self {
            Self::Abstract(s) => *s = "(-".to_owned() + s + ")",
            Self::Concrete(x) => *x *= -1.,
        }
    }
    pub fn simplify(&self) -> (Self, i32) {
        match self {
            Self::Abstract(_) => (self.clone(), 0),
            Self::Concrete(x) => {
                let y = x % (2. * std::f64::consts::PI);
                let k = y.div_euclid(std::f64::consts::PI / 2.);
                let rem = y % (std::f64::consts::PI / 2.);
                return (Self::Concrete(rem), k as i32);
            }
        }
    }
    pub fn is_zero_mod_two_pi(&self) -> bool {
        match self {
            Self::Abstract(_) => false,
            Self::Concrete(x) => (x % (2. * std::f64::consts::PI)).abs() < 1e-6,
        }
    }
    pub fn from_string(expr: String) -> Self {
        let as_f64 = expr.parse::<f64>();
        match as_f64 {
            Ok(x) => Self::Concrete(x),
            _ => Self::Abstract(expr),
        }
    }
    pub fn to_abstract(&self) -> Self {
        match self {
            Self::Abstract(x) => Self::Abstract(x.clone()),
            Self::Concrete(x) => Self::Abstract(x.to_string()),
        }
    }
    pub fn to_string(&self) -> String {
        match self {
            Self::Abstract(x) => x.clone(),
            Self::Concrete(x) => x.to_string(),
        }
    }
}
impl ops::Add<Parameter> for Parameter {
    type Output = Parameter;

    fn add(self, _rhs: Parameter) -> Parameter {
        match (&self, &_rhs) {
            (Self::Concrete(v1), Self::Concrete(v2)) => {
                Self::Concrete((v1 + v2) % (2. * std::f64::consts::PI))
            }
            _ => {
                let mut new_expr = self.to_string();
                new_expr.push_str("+");
                new_expr.push_str(&_rhs.to_string());
                Self::Abstract(new_expr)
            }
        }
    }
}
impl ops::AddAssign<Parameter> for Parameter {
    fn add_assign(&mut self, _rhs: Self) {
        match (&self, &_rhs) {
            (Self::Concrete(v1), Self::Concrete(v2)) => {
                *self = Self::Concrete((v1 + v2) % (2. * std::f64::consts::PI));
            }
            _ => {
                let mut new_expr = self.to_string();
                new_expr.push_str("+");
                new_expr.push_str(&_rhs.to_string());
                *self = Self::Abstract(new_expr);
            }
        }
    }
}
