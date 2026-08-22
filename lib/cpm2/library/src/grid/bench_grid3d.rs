fn bench_vector_push() {
    b.iter(|| {
        let mut vec = Vec::with_capacity(100);
        for i in 0..100 {
            vec.push(i);
        }
    });
}
