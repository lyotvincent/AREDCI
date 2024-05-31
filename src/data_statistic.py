# run in python3

# draw relation between loop length and #PET (supportive pair-end tags)
def draw_relation_between_loop_length_and_PET():
    import matplotlib.pyplot as plt
    import pandas as pd

    # read data
    data = open("D:/workspace/CHASOS/data/ChIA-PET2_result/GM12878_ChIA-PET2_result/GM12878.significant.interactions.intra.filtered.MICC", 'r')
    lines = data.readlines()
    data.close()
    length, pet = list(), list()
    for line in lines:
        line = line.strip().split('\t')
        length.append((int(line[4])+int(line[5])) - (int(line[1])+int(line[2])))
        pet.append(int(line[10]))
    # put length and pet into pandas dataframe
    # df = pd.DataFrame({'length': length, 'pet': pet})
    # draw
    plt.figure(figsize=(10, 5))
    plt.scatter(length, pet, s=10, c="r", marker="o", alpha=0.1)
    plt.xlabel("loop length")
    plt.ylabel("#PET")
    plt.xlim(0, 100000)
    plt.title("Relation between loop length and #PET")
    # plt.savefig("result/loop_length_and_PET.png")
    plt.show()


if __name__ == "__main__":
    draw_relation_between_loop_length_and_PET()

