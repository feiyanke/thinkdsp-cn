第九章：微分和积分
========================

这一章中，我们会继续分析时域中的窗函数和频域中的滤波的关系。
更具体的说，我们会分析有限差分（近似于微分）以及累加（近似于积分）两种操作的作用。

这章的代码 ``chap09.ipynb`` 可以在本书的 `代码库`_ 中找到，你也可以在 http://tinyurl.com/thinkdsp09 查看。

.. _代码库: https://github.com/AllenDowney/ThinkDSP

9.1 有限差分
---------------

在 :ref:`8.1 平滑` 中，我们对Fackbook的股价数据运用了一个平滑操作，发现在时域上平滑窗的作用
等价于在频域上的低通滤波。

在这一节中，我们来尝试计算它们每日价格差，你会发现在频域上，这相当于高通滤波的作用。

下面的代码从文件中读取数据保存到波形对象中并计算了它的频谱::

    import pandas as pd

    names = ['date', 'open', 'high', 'low', 'close', 'volume']
    df = pd.read_csv('fb.csv', header=0, names=names)
    ys = df.close.values[::-1]
    close = thinkdsp.Wave(ys, framerate=1)
    spectrum = wave.make_spectrum()

代码使用Pandas读取了数据，结果为 ``DataFrame`` 对象 ``df`` ，它的每一列代表了一天的开盘价，收盘价和最高最低价。
我选择了收盘价作为我们分析的数据，并把它保存到了波形对象中，这里的采样率为每天1次。

.. _图9.1:

.. figure:: images/thinkdsp048.png
    :alt: Daily closing price of Facebook and the spectrum of this time series
    :align: center

    图9.1： Fackbook股票收盘价的波形和频谱图

`图9.1`_ 展示了这些数据的波形和频谱图。直观上看，这个波形类似 :ref:`4.3 <4.3 布朗噪声>` 中的布朗噪声。
它的频谱看起来近似一条带有噪声的直线，斜率大约是-1.9， 与布朗噪声一样，它的斜率是常值。

接下来，我们使用 ``np.diff`` 来计算它的每日价格变化情况::

    diff = np.diff(ys)
    change = thinkdsp.Wave(diff, framerate=1)
    change_spectrum = change.make_spectrum()

`图9.2`_ 展示了结果的波形和频谱。可见，每日的变化曲线类似于白噪声，频谱近似一条平的直线，斜率大概是-0.06，接近于0，
这和白噪声是一样的。

.. _图9.2:

.. figure:: images/thinkdsp049.png
    :alt: Daily price change of Facebook and the spectrum of this time series
    :align: center

    图9.2： Fackbook股票每日价格变化的波形和频谱图

9.2 频域
--------------

实际上，计算相邻元素的差值与差分窗 *[1,-1]* 的卷积是一样的。你有可能会觉得这两个元素的值反了，实际上没有，
记住卷积的操作会先对窗函数进行反向的操作。

我们计算出了这个窗的DFT来看看它在频域上的作用效果::

    diff_window = np.array([1.0, -1.0])
    padded = thinkdsp.zero_pad(diff_window, len(close))
    diff_wave = thinkdsp.Wave(padded, framerate=close.framerate)
    diff_filter = diff_wave.make_spectrum()

如 `图9.3`_ 左图所示，有限差分窗对应了一个高通滤波器，它的幅值随着频率的增大而增大，
在低频的时候是线性的，而随着频率增大，曲线慢慢变平。下一节中，我们会分析这个原因。

.. _图9.3:

.. figure:: images/thinkdsp049.png
    :alt: Filters corresponding to the diff and differentiate operators (left) 
        and integration operator (right, log-y scale)
    :align: center

    图9.3： 差分和微分操作（左）以及积分操作（右）图的滤波效果

9.3 微分
-------------






























