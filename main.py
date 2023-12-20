from abstochkin.base import AbStochKin
import matplotlib.pyplot as plt

if __name__ == '__main__':
    sim = AbStochKin()
    sim.add_process_from_str("A -> B", 0.3, catalyst='E', Km=10)
    all_kwargs = [{"p0": {'A': i, 'B': 0, 'E': 10}, "t_max": 10} for i in range(40, 51)]
    sim.simulate_series_in_parallel(all_kwargs)

    plt.show()
