#generate plots
import matplotlib.pyplot as plt
import os

def plot_Losses(iter: list, d_loss: list, g_loss: list, mse_loss: list, title: str):

    #fig = plt.figure()
    plt.plot(iter, d_loss, "-b", label="d loss")
    plt.plot(iter, g_loss, "-r", label="g loss")
    plt.plot(iter, mse_loss, "-g", label="mse loss")

    plt.xlabel("n iteration")
    plt.legend(loc="upper left")
    plt.title(title)
    #plotFileName = os.path.join(pathSave, title + ".png")

    # save image
    #plt.savefig(plotFileName)  # should before show method

    # show
    plt.show()

    #return fig