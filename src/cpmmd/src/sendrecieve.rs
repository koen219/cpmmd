use serde_json::json;
use std::io::{Read, Write};
use std::net::TcpStream;

fn send_large_data(stream: &mut TcpStream, data: &[u8]) -> std::io::Result<()> {
    // Send the size of the data as a 4-byte header
    let size = data.len() as u32;
    stream.write_all(&size.to_be_bytes())?;
    // Send the actual data
    stream.write_all(data)?;
    Ok(())
}

pub fn send_json_message(
    stream: &mut TcpStream,
    message: &serde_json::Value,
) -> std::io::Result<()> {
    let data = message.to_string();
    send_large_data(stream, data.as_bytes())
}

pub fn receive_large_data(stream: &mut TcpStream) -> std::io::Result<Vec<u8>> {
    // Read the 4-byte size header
    let mut size_buf = [0; 4];
    stream.read_exact(&mut size_buf)?;
    let size = u32::from_be_bytes(size_buf) as usize;

    // Read the data in chunks
    let mut data = vec![0; size];
    stream.read_exact(&mut data)?;
    // let mut bytes_read = 0;
    //    while bytes_read < size {
    //        bytes_read += stream.read(&mut data[bytes_read..])?;
    //    }
    Ok(data)
}

pub fn receive_json_message(stream: &mut TcpStream) -> std::io::Result<serde_json::Value> {
    let data = receive_large_data(stream)?;
    let message: serde_json::Value = serde_json::from_slice(&data).expect("Failed to parse JSON");
    Ok(message)
}
