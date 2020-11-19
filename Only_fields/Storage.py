import pickle
import os

def save_data(Data, save_dir, filename):

    current_dir = os.getcwd()

    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)
        
    os.chdir(save_dir)
    pickle.dump( Data, open(filename + ".p", "wb" ) )
    os.chdir(current_dir)

def load_data(save_dir, filename):

    current_dir = os.getcwd()
    os.chdir(save_dir)
    load = pickle.load( open(filename + ".p", "rb" ) )
    os.chdir(current_dir)

    return load






















