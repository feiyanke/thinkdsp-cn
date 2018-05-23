第五章：自相关
===================

上一章中，我们说白噪声是不相关的，它任意时刻的值相互之间都是独立的。
而布朗噪声是相关的，它某个时刻的值依赖于上一个时刻。
这章中，我们引入 **自相关函数（autocorrelation function）** 
来精确的定义这些概念，这在信号分析中是很有用的一种方法。

这章的代码 ``chap05.ipynb`` 可以在本书的 `代码库`_ 中找到，你也可以在 http://tinyurl.com/thinkdsp05 查看。

.. _代码库: https://github.com/AllenDowney/ThinkDSP

5.1 相关性
---------------

通常来说，如果一个变量与另一个变量相关，那么意味着如果知道一个变量的值就可以得到另一个变量的一些信息。
我们有好几种方式来度量相关性，
其中最常见的一种是皮尔森相关系数（Pearson product-moment correlation coefficient），
记为： :math:`\rho` 。对于两个均包含N个测量值的变量 *x* 和 *y* ，

.. math::

    \rho  = \frac{{\sum\limits_i {({x_i} - {\mu _x})({y_i} - {\mu _y})} }}{{N{\sigma _x}{\sigma _y}}}

上式中的 :math:`{\mu _x}` 和 :math:`{\mu _y}` 是 *x* 和 *y* 的均值， :math:`{\sigma _x}` 和
 :math:`{\sigma _x}` 是它们的标准差。

皮尔森相关系数总是在[-1,1]的范围内，如果是正值，表示两个变量正相关，也就是说一个变量如果比较大，
那么另一个也趋向于变大，相反，如果是负值，那么表示两个变量如相关，也就是说一个变量如果比较大，
那么另一个会趋向于变小。

:math:`\rho` 的大小表示两个变量之间的相关程度，越大就越相关。如果 :math:`\rho` 是1或-1，
则说明它们是完全相关的，如果我们知道其中一个的值，就能准确预测出另一个值。如果 :math:`\rho` 
接近0，说明它们之间的相关性可能很弱，也就是知道其中一个的值，对于我们预测另一个的值来说没有什么意义。

这里，我用了“可能很弱”，是因为如果它们之间是非线性的关系，那么相关系数并不能很好的表征这种非线性关系。
在统计学中，非线性相关性也很重要，但是在信号处理中没有那么常见，因此这里我们先不考虑它。

Python中有很好几种方法来计算相关性。其中 ``np.corrcoef`` 可以计算多个变量两两之间的相关系数，结果用
一个 **相关矩阵（correlation matrix）** 来表示。

这里以两个变量为例，首先我构造了一个生成不同初始相位的正弦信号的函数，如下::

    def make_sine(offset):
        signal = thinkdsp.SinSignal(freq=440, offset=offset)
        wave = signal.make_wave(duration=0.5, framerate=10000)
        return wave

然后我们生成了两个不同相位的正弦信号::

    wave1 = make_sine(offset=0)
    wave2 = make_sine(offset=1)

`图5.1`_ 显示了这两个信号的前几个周期的波形。可以看到，一个信号的值比较大的时间，另一个的值也相对比较大，
感觉它们之间应该是相关的。

.. _图5.1:

.. figure:: images/thinkdsp026.png
    :alt: Two sine waves that differ by a phase offset of 1 radian; 
        their coefficient of correlation is 0.54
    :align: center

    图5.1： 两个相位差1弧度的正弦信号波形

它们之间的相关系数计算如下::

    >>> corr_matrix = np.corrcoef(wave1.ys, wave2.ys)
    [[ 1.    0.54]
    [ 0.54  1.  ]]

结果是一个相关矩阵，每个位置的值都代表的是对应的行号和列号两个信号之间的相关系数。
例如，第一行第一列元素的值代表了 ``wave1`` 和自身的相关系数，同样第二行第二列的值代表了
``wave2`` 和自身的相关系数。

其中分对角元素的值是我们所关心的，也就是 ``wave1`` 和 ``wave2`` 之间的相关系数0.54，
表明了它们之间有一定的相关性。

当它们之间的相位差增大时，它们的相关系数会减小直到相位差为180°，这时相关系数为-1。然后，
随着相位差继续增大，相关系数值又会增大直到360°，这个时候它们的相关系数为1.

`图5.2`_ 显示了它们之间的相关系数随相位差变化的情况。这个曲线的形状你应该很熟悉，它是余弦曲线。

.. _图5.2:

.. figure:: images/thinkdsp027.png
    :alt: The correlation of two sine waves as a function of 
        the phase offset between them. The result is a cosine
    :align: center

    图5.2： 正弦信号相关系数随相位差变化曲线

``thinkdsp`` 中提供一个简单的方法来计算波形之间的相关性::

    >>> wave1.corr(wave2)
    0.54

5.2 序列相关性
----------------

信号一般情况下都是对某个量的在一段时间内的连续测量值。
例如，声音信号就是在一段时间内我们对声压相对应的电压（电流）的测量值。

像这样的测量值都具有序列相关性，它表示一段测量序列中的不同时刻的测量值之间的相关性。
我们可以通过对信号进行移位后，在于原始信号计算相关系数来得到序列相关性::

    def serial_corr(wave, lag=1):
        n = len(wave)
        y1 = wave.ys[lag:]
        y2 = wave.ys[:n-lag]
        corr = np.corrcoef(y1, y2, ddof=0)[0, 1]
        return corr

``serial_corr`` 接受一个波形对象和 ``lag`` （一个整数，表示移位的位置）作为参数，
然后计算出移位后的波形和原始波形间的相关系数。

我们用上一章中的UG噪声信号来对这个函数进行测试::

    signal = thinkdsp.UncorrelatedGaussianNoise()
    wave = signal.make_wave(duration=0.5, framerate=11025)
    serial_corr(wave)

这个代码运行的结果是0.006，这表示信号的序列相关性很小。你可以试试用不同的 ``lag`` 值，
你会发现它们同样都很小。这是因为UG噪声是不相关的。

而对于布朗噪声来说，它的值是上一时刻的值加一个随机的步长，因此应该具有较大的序列相关性::

    signal = thinkdsp.BrownianNoise()
    wave = signal.make_wave(duration=0.5, framerate=11025)
    serial_corr(wave)

 我可以保证上面这段代码运行的结果要大于0.999.

 对于粉红噪声来说，它介于UG噪声和布朗噪声之间，相关系数也应该位于中间::

    signal = thinkdsp.PinkNoise(beta=1)
    wave = signal.make_wave(duration=0.5, framerate=11025)
    serial_corr(wave)

当 :math:`\beta = 1` 的时候，它的序列相关系数为0.851；
当 :math:`\beta = 0` 的时候，即UG噪声；
当 :math:`\beta = 2` 的时候，即布朗噪声。
如 `图5.3`_ 所示，粉红噪声的序列相关系数的范围为[0,1]。

.. _图5.3:

.. figure:: images/thinkdsp028.png
    :alt: Serial correlation for pink noise with a range of parameters
    :align: center

    图5.3： 粉红噪声的序列相关系数

5.3 自相关性
---------------

上一小节中我们计算了相邻两个时刻的信号值之间的相关性（也就是 ``lag=1`` ）。
实际上，我们可以很容易的计算出不同 ``lag`` 的相关系数。

实际上，你可以把 ``serial_corr`` 理解为一个从 ``lag`` 到相关系数的映射。
我们可以循环的计算出不同 ``lag`` 的相关系数::

    def autocorr(wave):
        lags = range(len(wave.ys)//2)
        corrs = [serial_corr(wave, lag) for lag in lags]
        return lags, corrs

``autocorr`` 接收一个波形对象作为参数，并返回一个序对形式的自相关函数：
``lags`` 是从0到波形长度一半的整数； ``corrs`` 是相对应的序列相关系数。

.. _图5.4:

.. figure:: images/thinkdsp029.png
    :alt: Autocorrelation functions for pink noise with a range of parameters
    :align: center

    图5.4： 不同参数的粉红噪声的自相关函数

`图5.4`_ 显示了三个不同 :math:`\beta` 的粉红噪声的自相关函数的曲线图。
可以看出，对于较小的 :math:`\beta` ，信号基本是不相关的，它的自相关函数很快的下降到0附近。
对于较大的 :math:`\beta` ，信号是强相关的，它的自相关函数比较大，并且下降比较慢。
当  :math:`\beta=1.7` 的时候，即使很大的 ``lag`` 值也有较强的相关性。这种情况，我们称为
**长期相关（long-range dependence）** ，因为信号的值依赖于长时间之前的值。

5.4 周期信号的自相关性
-----------------------





