import os

for dirname in os.listdir('processed_structures'):
    if os.path.isdir(os.path.join('processed_structures', dirname)):
        print('#'*20)
        print('Modeling ',dirname)
        print('#'*20)
        os.chdir('../STEGG_controler')
        os.system('python3 model_complex.py ../STEGG_benchmark_dataset/processed_structures/'+dirname+'/'+dirname+'.json')
        os.chdir('../STEGG_benchmark_dataset')
