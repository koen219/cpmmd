use std::process::Command;

fn main() {
    let status = Command::new("python3")
        .arg("gen_config.py")
        .status()
        .expect("Failed to run config generator");

    if !status.success() {
        panic!("Config generation failed");
    }

    // Make Cargo rebuild if these change
    println!("cargo:rerun-if-changed=src/config_spec.yaml");
    println!("cargo:rerun-if-changed=gen_config.py");
}
