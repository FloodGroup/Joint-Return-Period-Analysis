# 从这里开始

不要在这个文件夹里找一个叫 main 的平铺脚本。三个海湾各有一个主程序，就在发布目录里：

```text
gba/run_analysis.py
san_francisco/run_analysis.py
tokyo/run_analysis.py
```

打开其中任何一个，都能看到同一次分析的五次调用。五步的代码在 `jrpa/` 下，一步一个文件夹：

| 顺序 | 文件夹 | 这一步做什么 | 你要打开的文件 |
|--:|:--|:--|:--|
| 01 | `jrpa/step01_prepare/` | 读观测，做质控，换算水位基准 | `run.py` |
| 02 | `jrpa/step02_events/` | 合并降雨，挑选独立的降雨–水位事件 | `run.py` |
| 03 | `jrpa/step03_annual/` | 拟合年最大降雨 S1 和年最大水位 S2 | `run.py` |
| 04 | `jrpa/step04_copula/` | 拟合事件的边缘分布和 Copula | `run.py` |
| 05 | `jrpa/step05_export/` | 计算 Kendall 设计值，写出表格和报告 | `run.py` |

`jrpa/main.py` 只是把这五步按顺序串起来，供命令 `jrpa run` 和 `jrpa demo` 使用。平时运行上面的区域主程序即可。

`jrpa/common/` 放五步都会用到的设置、读写、基准和质控。它不是第六步。

研究区地图已经画好，在各区域的 `figures/` 里。画这些图的程序不在发布目录中。
