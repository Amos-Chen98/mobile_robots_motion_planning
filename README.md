# 深蓝学院《移动机器人运动规划》作业代码
课程地址：https://www.shenlanxueyuan.com/course/484

代码仓库：https://github.com/Amos-Chen98/mobile_robots_motion_planning

## hw1: Quick start

创建工作空间并拷贝相应功能包后，在`src/grid_path_searcher/launch/demo.launch`中，加一行

```xml
<node pkg="rviz" type="rviz" name="rivz" args="-d $(find grid_path_searcher)/launch/rviz_config/demo.rviz" />
```

然后运行此launch文件：

![](https://raw.githubusercontent.com/Amos-Chen98/Image_bed/main/2022/DeepinScreenshot_20220703113354.png)
