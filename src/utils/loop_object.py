

class Anchor():
    _instances = dict()

    def __new__(cls, start, end):
        key = (start, end)
        if key not in cls._instances:
            instance = super(Anchor, cls).__new__(cls)
            cls._instances[key] = instance
            instance.anchor_id = len(cls._instances)
            instance.start = start
            instance.end = end
        return cls._instances[key]

    def __init__(self, start, end):
        # self.anchor_name = f"anchor_{len(self._instances)}"
        # self.start = start
        # self.end = end
        pass

    def __str__(self):
        return f"({self.anchor_id},{self.start},{self.end})"

    __repr__ = __str__

    def __eq__(self, other):
        if isinstance(other, Anchor):
            return self.anchor_id == other.anchor_id and self.start == other.start and self.end == other.end
        return False

    def __hash__(self):
        return hash((self.anchor_id))


class Loop:
    _instances = dict()

    def __new__(cls, chr, anchor1, anchor2):
        key = (anchor1, anchor2)
        if key not in cls._instances:
            instance = super(Loop, cls).__new__(cls)
            cls._instances[key] = instance
            instance.loop_id = len(cls._instances)
            instance.chr = chr
            instance.anchor1 = anchor1
            instance.anchor2 = anchor2
        return cls._instances[key]

    def __init__(self, chr, anchor1, anchor2):
        # if (anchor1, anchor2) not in self._instances:
        #     self.loop_name = f"loop_{len(self._instances)}"
        # else:
        #     print("Warning: loop already exists, use the same loop name")
        #     print(self._instances[(anchor1, anchor2)])
        #     self.loop_name = self._instances[(anchor1, anchor2)].loop_name
        # self.anchor1 = anchor1
        # self.anchor2 = anchor2
        pass

    def __str__(self):
        return f"[{self.loop_id},{self.chr},{self.anchor1},{self.anchor2}]"

    __repr__ = __str__


    def __eq__(self, other):
        if isinstance(other, Loop):
            return self.loop_id == other.loop_id and self.chr == other.chr and self.anchor1 == other.anchor1 and self.anchor2 == other.anchor2
        return False

    def __hash__(self):
        return hash((self.loop_id))

    def get_loop_length(self):
        return abs((self.anchor2.start+self.anchor2.end)/2 - (self.anchor1.start+self.anchor1.end)/2)

    def overlap(self, other):
        if isinstance(other, Loop):
            return self.chr == other.chr and not (self.anchor1.start > other.anchor2.end and self.anchor2.start > other.anchor1.end)
        return False

class LoopInfo:
    def __init__(self, loop: Loop, pet, pvalue, fdr):
        self.loop = loop
        self.pet = pet
        self.pvalue = pvalue
        self.fdr = fdr

    def __str__(self):
        return f"<{self.loop},{self.pet},{self.fdr}>"

    __repr__ = __str__

    # we need to know if two loops are on the same anchors,
    # therefore, just compare loop, ignore chr, pet, fdr
    def __eq__(self, other):
        if isinstance(other, LoopInfo):
            return self.loop == other.loop
        elif isinstance(other, Loop):
            return self.loop == other
        return False

    def __hash__(self):
        return hash((self.loop.loop_id))

    def get_loop_id(self):
        return self.loop.loop_id

    def overlap(self, other):
        if isinstance(other, LoopInfo):
            return self.loop.overlap(other.loop)
        elif isinstance(other, Loop):
            return self.loop.overlap(other)
        return False

if __name__ == "__main__":
    print([Anchor(1, 2), 1, "tesstring"])
    print(Anchor(1, 2) == Anchor(1, 2))
    print(Anchor(1, 2) == Anchor(1, 3))
    anchors = [Anchor(1, 2), Anchor(1, 3), Anchor(2, 4)]
    import time
    time1 = time.time()
    print(Anchor(1, 2) in anchors)
    print(time.time() - time1)
    anchor1 = Anchor(1, 2)
    anchor2 = Anchor(1, 2)
    anchor3 = Anchor(1, 3)
    print("1", anchor1 == anchor2)
    # anchor1 and anchor2 is one object
    print("1", anchor1 is anchor2)
    print("2", anchor1 == anchor3)
    # anchor1 and anchor3 is not one object
    print("2", anchor1 is anchor3)
    print(anchor1.__hash__(), anchor2.__hash__(), anchor3.__hash__())
    loop1 = Loop(anchor1, anchor2)
    print("2.5", loop1)
    loop2 = Loop(anchor1, anchor2)
    loop3 = Loop(anchor1, anchor3)
    print("2.5", loop1)
    print("2.5", loop2)
    print("2.5", loop3)
    print("3", loop1 == loop2)
    # loop1 and loop2 is one object
    print("3", loop1 is loop2)
    print("4", loop1 == loop3)
    # loop1 and loop3 is not one object
    print("4", loop1 is loop3)
    print(loop1.__hash__(), loop2.__hash__(), loop3.__hash__())
    loop_info1 = LoopInfo(loop1, "chr1", 100, 0.01)
    loop_info2 = LoopInfo(loop2, "chr1", 80, 0.21)
    print("5", loop_info1 == loop_info2)
    loop_info_list = [loop_info1, loop_info2]
    loop_info3 = LoopInfo(loop3, "chr1", 233, 0.31)
    print("6", loop_info1 in loop_info_list)
    print("6", loop_info2 in loop_info_list)
    print("7", loop_info3 in loop_info_list)
    # print("11.5", loop_info_list.index(loop3))
    loop_info_list = [loop_info1, loop_info2, loop_info3]
    print("8", loop_info_list)
    print("9", set(loop_info_list))
    print("10", loop1 in loop_info_list)
    print("11", loop_info_list.index(loop3))
    loop_info4 = LoopInfo(loop1, "chr1", 100, 0.01)
    print("12", loop_info1)
    print("12", loop_info2)
    print("12", loop_info3)
    print("12", loop_info4)
    loop_info5 = LoopInfo(Loop(anchor1, anchor2), "chr1", 100, 0.01)
    print("13", loop_info5)
    loop_info6 = LoopInfo(Loop(anchor1, anchor2), "chr1", 100, 0.01)
    print("14", loop_info6)
    print(Loop._instances)
