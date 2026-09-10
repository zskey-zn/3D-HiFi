#!/usr/bin/awk -f

# 初始化输出列分隔符为制表符（也可以在命令行通过 -v OFS="\t" 传入）
BEGIN {
    OFS = "\t"
}

# 读取第一个文件：collapsed.contig.list
NR == FNR {
    # 存储每个 ID 对应的所有副本后缀
    for (i = 2; i <= NF; i++) {
        a[$1][i] = $i
    }
    next
}

# 读取第二个文件：merged_nodups.txt
NR > FNR {
    # 备份原始列，用于还原
    orig2 = $2
    orig6 = $6

    # 情况 1：$2 和 $6 都能在列表中匹配（执行两两组合）
    if ((orig2 in a) && (orig6 in a)) {
        for (i in a[orig2]) {
            for (j in a[orig6]) {
                $2 = orig2 "_" a[orig2][i]
                $6 = orig6 "_" a[orig6][j]
                print $0
            }
        }
    }

    # 情况 2：仅 $2 能在列表中匹配（$6 保持原样）
    if (orig2 in a) {
        for (i in a[orig2]) {
            $2 = orig2 "_" a[orig2][i]
            $6 = orig6
            print $0
        }
    }

    # 情况 3：仅 $6 能在列表中匹配（$2 保持原样）
    if (orig6 in a) {
        for (j in a[orig6]) {
            $2 = orig2
            $6 = orig6 "_" a[orig6][j]
            print $0
        }
    }

    # 情况 4：对应原本的 cat 逻辑，任何数据都无条件原样输出一份
    $2 = orig2
    $6 = orig6
    print $0
}
