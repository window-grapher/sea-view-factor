import numpy as np

def svf(terrain, obstacles, observer_pos, delta_angle, max_distance, print_info=False):
    """
    Calculates the Sea View Factor (SVF), which is the proportion of the sea visible from a given observer position.
    Parameters:
        terrain (np.ndarray): A 2D array representing the terrain elevation.
        obstacles (np.ndarray): A 2D array representing the height of obstacles on the terrain.
        observer_pos (tuple): A tuple (x, y) representing the observer's position on the terrain.
        delta_angle (float): The angular resolution in degrees for the visibility calculation.
        max_distance (int): The maximum distance to consider for visibility.
    CRS:
        EPSG: 4326
    Returns:
        tuple: A tuple containing:
            - sea_view_factor (float): The average proportion of the sea visible from the observer's position.
            - visibility (np.ndarray): An array of sea visibility fractions for each direction.
    """

    try:
        if len(terrain.shape) == 2:
            size_y, size_x = terrain.shape # 行数（y方向）, 列数（x方向）を別々に設定するように変更
        else:
            raise ValueError("[svf.py]Unexpected terrain shape: {}".format(terrain.shape))
        if print_info:
            print(f'[svf.py]size_x: {size_x}, size_y: {size_y}')

        observer_x, observer_y = observer_pos
        observer_height = terrain[observer_y, observer_x] + obstacles[observer_y, observer_x]
        if print_info:
            print(f'[svf.py]observer_x: {observer_x}, observer_y: {observer_y}, observer_height: {observer_height}')
        num_angles = int((2 * np.pi) / np.deg2rad(delta_angle))
        visibility = np.zeros(num_angles)

        # 観測者の視線の方向を0.5π（真北, 90度）スタートから-3π/2（真北, -270度）までの範囲で時計回りに設定
        angles = np.linspace(0.5 * np.pi, -1.5 * np.pi, num_angles, endpoint=False)

        for i, angle in enumerate(angles):
            step_x = np.cos(angle)
            step_y = np.sin(angle)
            if print_info:
                print(f'[svf.py]angle: {np.degrees(angle)}°, step_x: {step_x}, step_y: {step_y}')

            max_angle_of_elevation = - np.pi / 2  # 各角度での最大仰角（初期値はarctan2が取りえない-π/2）

            for distance in range(1, max_distance):
                ix = round(observer_x + step_x * distance)
                iy = round(observer_y - step_y * distance) # y座標は上下逆転しているのでプラスではなくマイナス

                if ix < 0 or ix >= size_x or iy < 0 or iy >= size_y: # グリッド外に出た場合
                    max_angle_of_elevation = 0  # 海に達する前にグリッド外に出たのでmax_angle_of_elevationの値は0
                    if print_info:
                        print(f'[svf.py]Out of Grid: {np.degrees(angle)}° at distance {distance} max_angle_of_elevation: {max_angle_of_elevation}')
                    break  # グリッド外に出た場合はその方向の探索を終了 # else文は実行されない

                height_at_point = terrain[iy, ix] + obstacles[iy, ix]
                relative_height = height_at_point - observer_height
                horizontal_distance = distance  # 水平距離としての距離

                if print_info:
                    print(f'[svf.py]distance: {distance}, ix: {ix}, iy: {iy}, height_at_point: {height_at_point}, relative_height: {relative_height}, horizontal_distance: {horizontal_distance}')

                angle_of_elevation = np.arctan2(relative_height, horizontal_distance) # angle_of_elevationはrelative_heightが負の場合は-π/2より大きく0の範囲
                if angle_of_elevation > max_angle_of_elevation:
                    max_angle_of_elevation = angle_of_elevation
                    if print_info:
                        print(f'[svf.py]max_angle_of_elevation updated: {max_angle_of_elevation}')

                if terrain[iy, ix] <= 0:  # max_angle_of_elevationは海に達した時までの最大値を保持 # ゼロメートル地帯の陸地の場合は要検討
                    if print_info:
                        print(f'[svf.py]reached Sea: {np.degrees(angle)}° at distance {distance} max_angle_of_elevation: {max_angle_of_elevation}')
                    # 海に達したらその方向の探索を終了 # else文は実行されない
                    # 手前に山がある地形で、海は見えないが奥に海がある場合、max_angle_of_elevationが0よりも大きい状況でBreakされる
                    break

            else: # for文が最後まで回りきった場合（海に達していない、かつグリッド外に出ていない）
                max_angle_of_elevation = 0 # 最終ループに達しても海ではなかったので、視野の範囲に海は入っていない
                if print_info:
                    print(f'[svf.py]No view of Sea at angle: {np.degrees(angle)}° max_angle_of_elevation: {max_angle_of_elevation}')

            if terrain[observer_y, observer_x] <= 0:  # ここまでで観測者が海にいる場合にはmax_angle_of_elevationが0になってしまう
                max_angle_of_elevation = - np.pi / 2  # 観測者が海にいる場合は天海率は1になるので代入しなおす
                if print_info:
                    print(f'[svf.py]Observer pos in Sea: {np.degrees(angle)}°')

            visible_AOE = max(0, - max_angle_of_elevation)
            # max_angle_of_elevation<0の場合は、海が見えるので、仰角（angle of elevation, AOE）の正負を逆転させて、天海率に変換。
            # max_angle_of_elevation>=0の場合は、海が見えないのでゼロとして扱う。

            # 天海率の計算（天空率から障害物による遮蔽を引く）
            sea_fraction = visible_AOE / (np.pi / 2)
            if print_info:
                print(f'[svf.py]sea_fraction: {sea_fraction} at {np.degrees(angle)}°')
            visibility[i] = sea_fraction  # -1より大きく1より小さい値

        # 全方向の平均を取って天海率を計算
        sea_view_factor = np.mean(visibility)
        return sea_view_factor, visibility

    except Exception as e:
        print(f"[svf.py]An error occurred during SVF calculations.: {e}")
        return None, None