# gen_config.py
import yaml
from pathlib import Path

TYPE_MAP_RUST = {
    'int': 'i32',
    'uint': 'usize',
    'float': 'f64',
    'List[float]': 'Vec<f64>',
    'bool': 'bool',
}

def gen_rust_struct(name, fields):
    out = "#[derive(Debug, Clone, serde::Deserialize)]\n"
    out += f"pub struct {name} {{\n"
    for key, t in fields.items():
        out += f"    pub {key}: {TYPE_MAP_RUST[t]},\n"
    out += "}\n"
    return out

def main():
    with open("config_spec.yaml") as f:
        spec = yaml.safe_load(f)

    for name, fields in spec.items():
        rust_code = gen_rust_struct(name, fields)
        Path("src").mkdir(exist_ok=True)
        with open(f"src/generated_config.rs", "w") as f:
            f.write(rust_code)

if __name__ == "__main__":
    main()
