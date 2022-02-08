#For R script

# Necessary packages
from modelV2 import Gain
from utils import plot_Losses
import utils
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

def main(args):
    '''Main function for UCI letter and spam datasets.

    Args:
      - data: .txt file
      - pathSave: directory to save plots
      - simulatedName: name of the plot figure
      - miss_rate: probability of missing components
      - batch:size: batch size
      - hint_rate: hint rate
      - alpha: hyperparameter
      - epochs: iterations

    Save:
    - imputed data
    - plot of the loss function
    '''

    pathSave = Path(args.pathSave)
    #data_path = Path(data_path)
    simulation = args.simulatedDataName
    #data_name = os.path.basename(data_path)
    data_name = args.data_name
    data = np.genfromtxt(data_name, delimiter=" ")
    #data = args.data

    gain_parameters = {'batch_size': args.batch_size,
                       'hint_rate': args.hint_rate,
                       'alpha': args.alpha,
                       'epochs': args.epochs}

    #Impute missing data
    imputed_data_x = Gain()
    imputed_data_x, d_loss, g_loss, mse_loss, it = imputed_data_x.runGain(data, gain_parameters)

    title = "Plot_"+simulation
    #generate plot
    plot_Losses(it, d_loss, g_loss, mse_loss, title, pathSave)

    #save numpy as csv
    #df = pd.DataFrame(imputed_data_x)
    #df.to_txt('imputedData.txt', index=False)

    #save numpy as txt
    np.savetxt("imputedData.txt", imputed_data_x, delimiter = " ")

    #return imputed_data_x


if __name__ == '__main__':
    # Inputs for the main function
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data_name",
        type=str)
    parser.add_argument(
        '--pathSave',
        type=str)
    parser.add_argument(
        "--simulatedDataName",
        type=str)
    parser.add_argument(
        '--miss_rate',
        help='missing data probability',
        default=0.2,
        type=float)
    parser.add_argument(
        '--batch_size',
        help='the number of samples in mini-batch',
        default=128,
        type=int)
    parser.add_argument(
        '--hint_rate',
        help='hint probability',
        default=0.9,
        type=float)
    parser.add_argument(
        '--alpha',
        help='hyperparameter',
        default=100,
        type=float)
    parser.add_argument(
        '--epochs',
        help='number of training interations',
        default=10000,
        type=int)

    args = parser.parse_args()

    # Calls main function
    main(args)

