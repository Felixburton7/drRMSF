import h5py

def inspect_h5_file(h5_file_path, domain_id, temperature, repl):
    with h5py.File(h5_file_path, 'r') as h5f:
        print(f"Available Temperatures: {list(h5f[domain_id].keys())}")
        if temperature not in h5f[domain_id]:
            print(f"Temperature {temperature}K not found.")
            return
        
        print(f"Available Replicas for {temperature}K: {list(h5f[domain_id][temperature].keys())}")
        for repl_id in h5f[domain_id][temperature].keys():
            print(f"\nReplica: {repl_id}")
            for dataset in h5f[domain_id][temperature][repl_id].keys():
                data = h5f[domain_id][temperature][repl_id][dataset]
                print(f"Dataset: {dataset}, Shape: {data.shape}, Type: {data.dtype}")
                for attr in data.attrs:
                    print(f"  Attribute - {attr}: {data.attrs[attr]}")
                
# Usage Example
inspect_h5_file(
    h5_file_path='mdcath_dataset_1r9lA02.h5',
    domain_id='1r9lA02',
    temperature='320',
    repl='0'
)
