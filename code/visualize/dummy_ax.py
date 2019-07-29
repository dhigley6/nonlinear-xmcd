"""Make dummy axes for labels/legends, etc.
"""

import matplotlib.pyplot as plt

def make_big_dummy_ax(fig):
    """Make a big dummy axis over the entire figure to use for global labels
    """
    dummyax = fig.add_subplot(111, frameon=False)
    dummyax.spines['top'].set_color('none')
    dummyax.spines['bottom'].set_color('none')
    dummyax.spines['left'].set_color('none')
    dummyax.spines['right'].set_color('none')
    dummyax.tick_params(labelcolor='none', top='off', left='off', right='off', bottom='off')
    #plt.setp(dummyax.get_xticklabels(), alpha=0)
    plt.setp(dummyax.get_yticklabels(), alpha=0)
    return dummyax
