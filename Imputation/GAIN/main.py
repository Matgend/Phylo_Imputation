from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from data_loader import data_loader
# Necessary packages
from model import Gain
import argparse
import numpy as np




def main(args):
    '''Main function for UCI letter and spam datasets.

    Args:
      - data_name: letter or spam
      - miss_rate: probability of missing components
      - batch:size: batch size
      - hint_rate: hint rate
      - alpha: hyperparameter
      - iterations: iterations

    Returns:
      - imputed_data_x: imputed data
      - rmse: Root Mean Squared Error
    '''

    data_name = args.data_name
    miss_rate = args.miss_rate

    """
    gain_parameters = {'batch_size': args.batch_size,
                       'hint_rate': args.hint_rate,
                       'alpha': args.alpha,
                       'iterations': args.iterations}
    """

    # Load data and introduce missingness
    ori_data_x, miss_data_x, data_m = data_loader(data_name, miss_rate)

    # Impute missing data
    imputed_data_x = Gain()
    imputed_data_x = imputed_data_x.runGain(ori_data_x, args.batch_size, args.hint_rate, args.alpha, args.iteration)
    np.savetxt("imputedData", imputed_data_x, delimiter=",")

    # Report the RMSE performance
    rmse = imputed_data_x._rmse_loss(ori_data_x, imputed_data_x, data_m)

    print()
    print('RMSE Performance: ' + str(np.round(rmse, 4)))

    return imputed_data_x, rmse


if __name__ == '__main__':
    # Inputs for the main function
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--data_name',
        choices=['letter', 'spam'],
        default='spam',
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
        '--iterations',
        help='number of training interations',
        default=10000,
        type=int)

    args = parser.parse_args()

    # Calls main function
    imputed_data, rmse = main(args)
