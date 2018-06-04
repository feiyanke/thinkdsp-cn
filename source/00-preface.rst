前言
========

信号处理是我很感兴趣的一个学科，它被广泛的应用在很多领域。
如果理解了信号处理中的一些基本概念，可以帮助我们更好的理解我们生活中所看到
和所听到的东西。

但是，除非你是学机械或电子专业的，你大概很难有机会去了解信号处理这个学科。
大多数的相关的书籍（包括学校使用的教材）都使用的是自底向上的形式，
先从最基本的数学抽象开始讲起（例如向量）然后逐渐深入，并且通常缺少具体的应用示例以及
与现实世界的对应关系。

阅读本书的前提是你需要有些编程的基础，这样你就可以可以使用这项技能来学习一些有趣的东西了。
使用编程式的讲解方式，我们就可以很直观的展示大多数的概念，比如在第一章结束后，你就可以通过程序来
分析声音以及其他的信号，生成新的声音了。每一章中，我们都会介绍一项新的技术，以及它在处理真实
信号中的应用。我们会先学习怎样使用这个技术，然后再学习它的工作原理。

我相信这样的学习形式是更实用也更有趣的（希望你们是这样认为的）。

本书的读者
-------------

本书使用Python作为编程语言，因此需要读者有一些Python的编程基础并且对面向对象的程序设计有所了解。

如果你对Python还不太熟的话，你可以先看看我写的另一本书，ThinkPython，它是一本Python的入门数据，
适用于那些没有编程基础的人。如果你有编程经验但是对Python不熟悉，那么你也可以看Mark Lutz的Learning Python。

本书中还会用到Numpy和Scipy两个Python的扩展库，如果你对此不熟悉也没有关系，我会在使用它们的函数和类的时候，
进行简单的解释。

本书的代码
--------------

本书的代码和书中使用到的声音文件都托管在Github上： https://github.com/AllenDowney/ThinkDSP 。
Git是一个版本管理工具，使用它你可以很好的管理和跟踪整个项目中的所有文件，被Git管理项目我们称之为代码库。
Github是一个托管代码库的网站，它提供了很方便的用户界面来管理你的代码库。

本书代码库的Github主页上提供了几种方式来使用其中的代码：

* 你可以把我的代码库Fork到自己的账号下。（如果你没有Github账户，你可以免费注册一个）。
  Fork后你就有了一个自己的代码库了，你可以使用它来跟踪管理你在学习过程中编写的代码。
  然后，你需要使用clone命令来把代码库中的代码复制到你的本地电脑中。

* 如果你没有Github账号，你也可以直接clone我的代码库。但是这样你就不能把你写的代码上传到Github上了，不过你
  依然可以在本地跟踪和管理你的代码。

* 如果你不想使用Git，你还可以把代码库下载成一个zip文件。（Github页面的右下角有一个下载的按钮）

本书中的所有代码均可以工作与Python 2和Python 3。

我在开发这些代码的时候，使用了Continuum Analytics公司开发的Anaconda，它是一个免费的Python发行版，其中包含了
我们将要用到的大部分Python包。Anaconda的安装非常简单，默认情况下它使用用户级的安装方式，因此不需要管理员权限。
它同时支持Python2和Python3。你可以在 http://continuum.io/downloads 下载Anaconda。

如果你不想使用Anaconda，你需要自己手动安装以下的Python包：

* Numpy：用于数值计算， http://www.numpy.org/ ；

* Scipy： 用于科学计算， http://www.scipy.org/ ；

* matplotlib：用于作图， http://matplotlib.org/ 。

虽然这些程序包都是很常用的，但是它们不会包含在Python的安装包中，因此你需要手动安装。
在某些环境下，它们有可能不太好安装，如果你在安装的时候遇到问题，我建议你使用Anaconda或其他的Python发行版，
一般这些发行版的Python都默认包含了这些程序包。

本书的练习有的需要使用Jupyter notebook，如果你对此不熟悉，可以参考 http://jupyter.org 。
使用Jupyter有三种方式：

**在电脑中运行Jupyter**

如果安装了Anaconda，其中默认包含了Jupyter，你可以使用命令行来启动Jupyter服务::

    $ jupyter notebook

如果没有安装Jupyter，可以使用下面的命令行进行安装::

    $ conda install jupyter

当你启动了Jupyter服务后，它会在浏览器中自动的打开一个新的页面也就是Jupyter notebook的主页。

**在Binder中运行Jupyter**

Binder是一个运行Jupyter的网络服务，通过 http://mybinder.org/repo/AllenDowney/ThinkDSP 你可以直接打开
本书的Jupyter主页，并且你可以在编辑和运行里面的代码，但是你编写和改动的代码不会被保存起来，如果你把页面关闭了或者
长时间没有操作（1个小时），那么这些代码会消失。

**在nbviewer中查看**

本书还会提供在noviewer中的链接，在nbviewer中你仅仅可以静态展示代码和结果。你可以通过nbviewer的链接来阅读本书的代码，
也可以播放示例中的声音，但是你不能进行改动也不能运行它们，交互式的控件也使用不了。

Good Luck，and have fun!

Contributor List
----------------

If you have a suggestion or correction, please send email to downey@allendowney.com. 
If I make a change based on your feedback, I will add you to the contributor list (unless you ask to be omitted).

If you include at least part of the sentence the error appears in, 
that makes it easy for me to search. Page and section numbers are fine, too, but not as easy to work with. Thanks!

* Before I started writing, my thoughts about this book benefited from conversations 
  with Boulos Harb at Google and Aurelio Ramos, formerly at Harmonix Music Systems.

* During the Fall 2013 semester, Nathan Lintz and Ian Daniher worked with me on an independent 
  study project and helped me with the first draft of this book.

* On Reddit’s DSP forum, the anonymous user RamjetSoundwave helped me fix a problem with 
  my implementation of Brownian Noise. And andodli found a typo.

* In Spring 2015 I had the pleasure of teaching this material along with Prof. 
  Oscar Mur-Miranda and Prof. Siddhartan Govindasamy. Both made many suggestions and corrections.

* Giuseppe Masetti sent a number of very helpful suggestions.

* Kim Cofer copyedited the entire book, found many errors, and made many helpful suggestions.

Other people who found typos and errors include Silas Gyger and Abe Raher.

Special thanks to the technical reviewers, Eric Peters, Bruce Levens, and John Vincent, 
for many helpful suggestions, clarifications, and corrections.

Also thanks to Freesound, which is the source of many of the sound samples I use in this book, 
and to the Freesound users who contributed those samples. 
I include some of their wave files in the GitHub repository for this book, using the original file names, 
so it should be easy to find their sources.

Unfortunately, most Freesound users don’t make their real names available, 
so I can only thank them by their user names. 
Samples used in this book were contributed by Freesound users: 
iluppai, wcfl10, thirsk, docquesting, kleeb, landup, zippi1, themusicalnomad, 
bcjordan, rockwehrmann, marcgascon7, jcveliz. Thank you all!