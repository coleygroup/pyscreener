from pathlib import Path
import subprocess as sp
import time

def main():
    for config in Path('test_configs').iterdir():
        print(sp.run([
            'sbatch',
            '-o', f'{config.stem}_%j.out',
            '-e', f'{config.stem}_%j.err', 
            'run_pyscreener.batch', str(config)]))
        time.sleep(1)

if __name__ == "__main__":
    main()