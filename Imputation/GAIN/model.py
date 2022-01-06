import numpy as np
from typing import List
from tqdm import tqdm
import tensorflow as tf
from sklearn.preprocessing import MinMaxScaler


class Gain:

    def __init__(self):
        self.scaler = MinMaxScaler()

    def normalization(self, data: np.array) -> np.array:
        """
        Arg:
            data: original data (numpy array)

        Output:
            normalized_data: numpy array
        """
        self.scaler.fit(data)
        normalized_data = self.scaler.transform(data)
        return normalized_data

    def denormalization(self, normalized_data: np.array) -> np.array:

        denormalized_data = self.scaler.inverse_transform(normalized_data)
        return denormalized_data

    def _xavier_init(self, shape: List[int]) -> tf.tensor:
        """Xavier initialization.

          Args:
            - size: vector size

          Returns:
            - initialized random vector.
          """
        in_dim = shape[0]
        xavier_stddev = 1. / tf.sqrt(in_dim / 2.)
        return tf.random_normal(shape=shape, stddev=xavier_stddev)

    def sample_batch_index(self, total: int, batch_size: int) -> int:
        """
        Sample index of the mini-batch.

        Args:
          - total: total number of samples
          - batch_size: batch size

        Returns:
          - batch_idx: batch index
        """
        total_idx = np.random.permutation(total)
        batch_idx = total_idx[:batch_size]
        return batch_idx

    def _rmse_loss(self, ori_data: np.array, imputed_data: np.array, data_m) -> float:
        """Compute RMSE loss between ori_data and imputed_data

        Args:
        - ori_data: original data without missing values
        - imputed_data: imputed data
        - data_m: indicator matrix for missingness (np array?

        Returns:
        - rmse: Root Mean Squared Error
        """
        ori_data = self.normalization(ori_data)
        imputed_data = self.denormalization(imputed_data)

        # Only for missing values
        nominator = np.sum(((1-data_m) * ori_data - (1-data_m) * imputed_data)**2)
        denominator = np.sum(1-data_m)

        rmse = np.sqrt(nominator/float(denominator))

        return rmse

    def train(self, data: np.array, batch_size: int, hint_rate: float, alpha: int, epochs: int) -> np.array:

        # normalize data
        normalized_data = self.normalization(data)
        # replace NaN by 0 in normalized data
        normalized_data_NaN0 = np.nan_to_num(normalized_data, 0)

        row, col = data.shape
        # GAIN architecture #
        # Input placeholders
        # Data vector
        X = tf.placeholder(tf.float32, shape=[None, col])
        # Mask vector
        M = tf.placeholder(tf.float32, shape=[None, col])
        # Hint vector
        H = tf.placeholder(tf.float32, shape=[None, col])

        # Discriminator variables
        D_W1 = tf.Variable(self._xavier_init([col * 2, col]))# Data + Hint as inputs
        D_b1 = tf.Variable(tf.zeros(shape=[col]))

        D_W2 = tf.Variable(self._xavier_init([col, col]))
        D_b2 = tf.Variable(tf.zeros(shape=[col]))

        D_W3 = tf.Variable(self._xavier_init([col, col]))
        D_b3 = tf.Variable(tf.zeros(shape=[col]))  # Multi-variate outputs

        theta_D = [D_W1, D_W2, D_W3, D_b1, D_b2, D_b3]

        # Generator variables
        # Data + Mask as inputs (Random noise is in missing components)
        G_W1 = tf.Variable(self._xavier_init([col * 2, col]))
        G_b1 = tf.Variable(tf.zeros(shape=[col]))

        G_W2 = tf.Variable(self._xavier_init([col, col]))
        G_b2 = tf.Variable(tf.zeros(shape=[col]))

        G_W3 = tf.Variable(self._xavier_init([col, col]))
        G_b3 = tf.Variable(tf.zeros(shape=[col]))

        theta_G = [G_W1, G_W2, G_W3, G_b1, G_b2, G_b3]

        ## GAIN functions
        # Generator
        def generator(x, m):
            # Concatenate Mask and Data
            inputs = tf.concat(values=[x, m], axis=1)
            # an entire layer
            G_h1 = tf.nn.relu(tf.matmul(inputs, G_W1) + G_b1)
            # second entire layer
            G_h2 = tf.nn.relu(tf.matmul(G_h1, G_W2) + G_b2)
            # MinMax normalized output
            G_prob = tf.nn.sigmoid(tf.matmul(G_h2, G_W3) + G_b3)
            return G_prob

        # Discriminator
        def discriminator(x, h):
            # Concatenate Data and Hint
            inputs = tf.concat(values=[x, h], axis=1)
            D_h1 = tf.nn.relu(tf.matmul(inputs, D_W1) + D_b1)
            D_h2 = tf.nn.relu(tf.matmul(D_h1, D_W2) + D_b2)
            D_logit = tf.matmul(D_h2, D_W3) + D_b3
            D_prob = tf.nn.sigmoid(D_logit)
            return D_prob

        ## GAIN structure
        # Generator
        G_sample = generator(X, M)

        # Combine with observed data
        Hat_X = X * M + G_sample * (1 - M)

        # Discriminator
        D_prob = discriminator(Hat_X, H)

        ## GAIN loss
        D_loss_temp = -tf.reduce_mean(M * tf.log(D_prob + 1e-8) \
                                      + (1 - M) * tf.log(1. - D_prob + 1e-8))

        G_loss_temp = -tf.reduce_mean((1 - M) * tf.log(D_prob + 1e-8))

        MSE_loss = \
            tf.reduce_mean((M * X - M * G_sample) ** 2) / tf.reduce_mean(M)

        D_loss = D_loss_temp
        G_loss = G_loss_temp + alpha * MSE_loss # learning rate

        ## GAIN solver
        D_solver = tf.train.AdamOptimizer().minimize(D_loss, var_list=theta_D)
        G_solver = tf.train.AdamOptimizer().minimize(G_loss, var_list=theta_G)

        ## Iterations
        sess = tf.Session()
        sess.run(tf.global_variables_initializer())

        # Start Iterations
        for it in tqdm(range(epochs)):
            # Sample batch
            batch_idx = sample_batch_index(row, batch_size)
            X_mb = normalized_data_NaN0[batch_idx, :] # mb = mini batch
            M_mb = data_m[batch_idx, :]
            # Sample random vectors
            Z_mb = uniform_sampler(0, 0.01, batch_size, col)
            # Sample hint vectors
            H_mb_temp = binary_sampler(hint_rate, batch_size, col)
            H_mb = M_mb * H_mb_temp

            # Combine random vectors with observed vectors
            X_mb = M_mb * X_mb + (1 - M_mb) * Z_mb

            _, D_loss_curr = sess.run([D_solver, D_loss_temp],
                                      feed_dict={M: M_mb, X: X_mb, H: H_mb})
            _, G_loss_curr, MSE_loss_curr = \
                sess.run([G_solver, G_loss_temp, MSE_loss],
                         feed_dict={X: X_mb, M: M_mb, H: H_mb})

        ## Return imputed data
        Z_mb = uniform_sampler(0, 0.01, row, col)
        M_mb = data_m
        X_mb = normalized_data_NaN0
        X_mb = M_mb * X_mb + (1 - M_mb) * Z_mb

        imputed_data = sess.run([G_sample], feed_dict={X: X_mb, M: M_mb})[0]

        imputed_data = data_m * normalized_data_NaN0 + (1 - data_m) * imputed_data

        # Renormalization
        imputed_data = renormalization(imputed_data, norm_parameters)

        # Rounding
        imputed_data = rounding(imputed_data, data_x)


    def test(self):
        pass

    def predict(self):
        pass

    def export(self):
        pass

    def load(self):
        pass

