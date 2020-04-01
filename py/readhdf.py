import h5py

with h5py.File('output_group.h5', 'r') as f:
    print(f['hoge/x'].value)  # => 200
    print(f['hoge/a'].value)

print(f)
