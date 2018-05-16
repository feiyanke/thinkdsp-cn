第三章：非周期信号
====================

到目前为止，我们学习到的信号都是周期的，它们会无限的重复并循环下去。
这也意味着，它们的频谱不会随着时间的变化而变化。
本章我们会学习非周期的信号，它们的频谱是随时间而变化。
几乎现实中的所有声音信号都是这样的非周期信号。

本章还会介绍声谱图（spectrograms），用于非周期信号的频谱的可视化。

这章的代码 ``chap03.ipynb`` 可以在本书的 `代码库`_ 中找到，你也可以在 http://tinyurl.com/thinkdsp03 查看。

.. _代码库: https://github.com/AllenDowney/ThinkDSP

3.1 线性啁啾声

我们先从啁啾声（ **chirp** ）开始，这是一种频率随时间变化的信号。 ``thinkdsp`` 中提供了一个 ``Chirp`` 类，
用来生成频率在一定范围内线性变化的正弦信号。

以下的代码生成了一个频率从220Hz变化到880Hz（从A3到A5升了两个八度）的 ``Chirp`` 信号::

    signal = thinkdsp.Chirp(start=220, end=880)
    wave = signal.make_wave()

`图3.1`_ 的三个图分别显示了这个信号开始，中间和结束的三段波形图。图中可以清晰的看出频率的变化。

.. _图3.1:

.. figure:: images/thinkdsp012.png
    :alt: Chirp waveform near the beginning, middle, and end
    :align: center

    图3.1： Chirp信号开始，中间和结束的波形图

我们先来看看 ``Chirp`` 是怎么实现的，下面是这个类的代码::

    class Chirp(Signal):
    
        def __init__(self, start=440, end=880, amp=1.0):
            self.start = start
            self.end = end
            self.amp = amp

        def evaluate(self, ts):
            freqs = np.linspace(self.start, self.end, len(ts)-1)
            return self._evaluate(ts, freqs)

        def _evaluate(self, ts, freqs):
            dts = np.diff(ts)
            dphis = PI2 * freqs * dts
            phases = np.cumsum(dphis)
            phases = np.insert(phases, 0, 0)
            ys = self.amp * np.cos(phases)
            return ys

构造函数 ``__init__`` 中 ``start`` 和 ``end`` 分别表示开始和结束的频率（单位为Hz）。 ``amp`` 表示幅值。

``evaluate`` 方法实现了信号的计算，其中 ``ts`` 表示采样点的时间序列。简单起见，我们假设采样的时间间隔是固定的。

假设 ``ts`` 的长度为 n, 就有 n-1 个时间段， 我们可以使用 ``np.linspace`` 求出对应的频率，也就是从 ``start``
到 ``end`` 的等间隔和 n-1 个值（Numpy数组）

``_evaluate`` 私有方法完成了接下来的数学运算  [1]_ 。其中 ``np.diff`` 计算了 ``ts`` 中相邻的采样点的时间间隔，
结果保存在 ``dfs`` 中。如果 ``ts`` 中的元素是等距的，那么 ``dts`` 中的值应该都是一样的。

接下来我们来计算在每个时间间隔内的相位变化。在 `1.7 信号对象`_ 中，我们看到当频率是常量的时候，相位 :math:`\varphi`
是相对时间线性变化的：

.. math::

    \varphi  = 2\pi ft

当频率随时间变化的时候，相位在短时间 :math:`\Delta t` 内的变化量是：

.. math::

    \Delta \varphi  = 2\pi f(t)\Delta t

因为 ``freqs`` 实际上就是 :math:`f(t)` ，而 ``dts`` 就是 :math:`\Delta t` ，因此上式可以写成如下的Python代码::

    dphis = PI2 * freqs * dts

现在 ``dphis`` 中包含了相位的变化量，我们可以通过累加来得到各个时间点的相位::

    phases = np.cumsum(dphis)
    phases = np.insert(phases, 0, 0)

``np.cumsum`` 方法计算出了累加值，可以看出来这个值的第一个元素不为0，因此我们需要使用 ``np.insert`` 在前面添加一个0值。
最后， 我们使用 ``np.cos`` 计算出了整个信号的值（记住相位是用弧度表示的）。

实际上，如果用微积分来表示，当 :math:`\Delta t` 足够小的时候：

.. math::

    d\varphi  = 2\pi f(t)dt

两边同时除以 :math:`dt` 得到：

.. math::

    \frac{{d\varphi }}{{dt}} = 2\pi f(t)

也就是说，相位的微分就是频率。反过来，相位应该是频率的积分。所以，我们可以使用 ``np.cumsum`` 累加来得到相位，
因为累加实际上就是积分的近似计算方法。

.. admonition:: 译者注

    译者觉得这里其实没必要使用累加来计算，因为线性变化的频率很容易可以求得解析解。考虑到相位的微分就是频率，
    如果希望频率线性变化，那么 :math:`f(t) = at + b`，
    那么相位就是频率的不定积分 :math:`\varphi  = \frac{a}{2}{t^2} + bt + c` ，
    考虑到单位转换后，最后的信号随时间变化公式为：:math:`s = \cos (2\pi (\frac{a}{2}t + b)t + c)`
    已知 ``ts`` 和 ``start`` ``end`` ，很容易可以得到 a 和 b 的值，最终 ``evaluate`` 的代码如下::

        def evaluate(self, ts):
            k = (self.end - self.start) / (ts[-1] - ts[0]) / 2
            freqs = k*(ts-ts[0])+self.start
            phases = 2 * np.pi * freqs * ts
            ys = self.amp * np.cos(phases)
            return ys

    有兴趣的读者可以尝试使用这个方式来构造信号，看和书中提供的方法产生的信号是否一致。

3.2 指数啁啾声
-----------------

当你听这个啁啾声的时候，你会发现一开始音高上升的很快，然后会慢下来。这个啁啾声跨越了两个八度，
跨越第一个八度只用了 1/3s 时间，而第二个八度用了 2/3s。

造成这个现象的原因是我们感受到的音高取决于频率的对数，也就是说我们听到的两个声音的音高间隔
取决于它们之间的频率比值，而不是差值。用音乐的术语来说，两个音高之间的间隔，被称为音程（ **interval** ）

例如，一个八度指的是频率之比为2的两个音高之间的间隔。因此从220Hz到440Hz为一个八度，
从440Hz到880Hz又是一个八度。虽然他们之间的频率差更大，但是他们的音程是一样的。

因此，如果频率是线性升高的，那么听起来音高是按对数升高的。

如果我们想得到音高按线性变化的信号，那么信号的频率就得按指数变化。这种信号我们成为指数啁啾声。代码如下::

    class ExpoChirp(Chirp):
    
        def evaluate(self, ts):
            start, end = np.log10(self.start), np.log10(self.end)
            freqs = np.logspace(start, end, len(ts)-1)
            return self._evaluate(ts, freqs)

这里我们使用了 ``np.logspace`` 来替代 ``np.linspace`` ，它可以产生按指数变化的序列值。

其他的代码与之前的 ``Chirp`` 是一样的，我们使用它来生成一个指数啁啾声::

    signal = thinkdsp.ExpoChirp(start=220, end=880)
    wave = signal.make_wave(duration=1)

你可以在 ``chap03.ipynb`` 中听一听这些信号的区别。

3.3 啁啾声的频谱
---------------------

啁啾声的频谱图是怎样的呢？这里我们构造了一个1s内八度的信号，并且计算出了它的频谱::

    signal = thinkdsp.Chirp(start=220, end=440)
    wave = signal.make_wave(duration=1)
    spectrum = wave.make_spectrum()

`图3.2`_ 展示了这个频谱图。可以看到，这个信号包含从220Hz到440Hz的所有频率成分。
还可以注意到，在220Hz到440Hz区间内，频谱图大概是平的，这就表明频率在时间上是均匀变化的。
基于此我们可以猜测指数啁啾声的频谱是什么样子的吗？

.. _图3.2:

.. figure:: images/thinkdsp013.png
    :alt: Spectrum of a one-second one-octave chirp
    :align: center

    图3.2： 1s升八度的Chirp信号

实际上，从频谱图中，我们可以得到信号的频率成分的信息，但是却掩盖了频率随时间变化的信息。
我们不能从频谱中看出信号的频率是随时间变大还是变小了。

3.4 声谱图
----------------

为了展示信号频率随时间变化的关系，我们可以把信号分段后分别计算频谱，然后画出每段的频谱图。
这种方法我们成为 **短时傅立叶变换（STFT）** 。

我们常用声谱图（ **spectrogram** ）来可视化STFT的结果。声谱图的x轴是时间，y轴是频率。
声谱图中的每列显示了一小段时间内信号的频谱，使用灰度值（或颜色亮度）来表示幅值大小。

我们以 ``Chirp`` 信号作为例子来计算声谱图::

    signal = thinkdsp.Chirp(start=220, end=440)
    wave = signal.make_wave(duration=1, framerate=11025)
    spectrogram = wave.make_spectrogram(seg_length=512)
    spectrogram.plot(high=700)

``Wave`` 类提供了 ``make_spectrogram`` 来生成声谱图。其中 ``seg_length`` 表示每段包含的采样点数。
这里使用了512，通常情况使用2的n次方的值可以提升FFT的效率。`图3.3`_ 为生成的声谱图。

.. _图3.3:

.. figure:: images/thinkdsp014.png
    :alt: Spectrogram of a one-second one-octave chirp
    :align: center

    图3.3： 1s升八度的Chirp信号的声谱图

图中，x轴的时间范围从0s到1s，y轴频率范围从0Hz到700Hz。因为信号频率成分比较低，为了更清除的展示，
我把整个声谱图的上部分裁剪了，实际上完整的频率范围是0~5512.5Hz，即采样率的一半。

声谱图清楚的展示了信号频率随时间的变化情况。但是，我们也可以注意到，图中每列的峰值都有2-3个单位的模糊，
实际上这反应了声谱图的频率分辨率是有限的。



    



.. [1] 方法名前面加下划线表示这个方法是私有的，不应该在外部进行调用
    


