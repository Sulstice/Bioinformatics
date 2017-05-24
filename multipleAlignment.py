import time
from sys import getsizeof
import math as math
import numpy as np

global MATCH
global MISMATCH
global INDEL
global scoredArray
global seqLen1
global seqLen2
global seqLen3
global maxValue

def createScoring():
    i = 0
    j = 0
    k = 0

    scoringArray = np.zeros((5, 5, 5))

    for i in range (0, 5):
        for j in range(0, 5):
            for k in range(0, 5):
                score = 0
                if (i == 4) and (j == 4) and (k == 4):
                    score = 0
                elif (i == 4 and j == 4 and k != 4) or (i == 4 and j != 4 and k == 4) or (i != 4 and j == 4 and k == 4):
                    score = INDEL * 2
                elif (i == 4 and j != 4 and k != 4):
                    if (j == k):
                        score = (INDEL * 2) + MATCH
                    else:
                        score = (INDEL * 2) + MISMATCH
                elif (i != 4 and j == 4 and k != 4):
                    if (i == k):
                        score = (INDEL * 2) + MATCH
                    else:
                        score = (INDEL * 2) + MISMATCH
                elif (i != 4 and j != 4 and k== 4):
                    if (i == j):
                        score = (INDEL * 2) + MATCH
                    else:
                        score = (INDEL * 2) + MISMATCH
                elif ( i < 4 and j < 4 and k < 4):
                    if (i == j and i == k):
                        score = MATCH * 3
                    elif ( i == j and i != k) or (i != j and i == k) or (j == k and j != i) or (j != k and j == i):
                        score = MATCH + (2 * MISMATCH)
                    elif (i != j and j != k and k != i):
                        score = MISMATCH * 3

                scoringArray[i, j, k] = score

    return scoringArray


def alignment(iStart, jStart, kStart, LengthSeq1, LengthSeq2, LengthSeq3):

    # Score Forward

    lengthi = LengthSeq1 - iStart
    lengthj = LengthSeq2 - jStart
    lengthk = LengthSeq3 - kStart

    print (lengthi)
    if (lengthi <= 1):
        return

    mid = math.floor((LengthSeq1 + iStart) / 2)

    score1 = np.zeros((2, int(lengthj + 1), int(lengthk + 1)))
    score2 = np.zeros((2, int(lengthj + 1), int(lengthk + 1)))

    for i in range(0, 2):
        for j in range(0, int(lengthj + 1)):
            for k in range(0, int(lengthk + 1)):
                score1[i, j, k] = -100000

    score1[0, 0, 0] = 0

    for k in range(1, lengthk + 1):
        score1[0, 0, k] = scoredArray[4, 4, seqLen3[kStart + k]] * k
    for j in range(1, lengthj + 1):
        score1[0, j, 0] = scoredArray[4, seqLen2[jStart + j], 4]

    for j in range(0, lengthj + 1):
        for k in range(0, lengthk + 1):
            scoreXYplane = 0
            subscore = np.array([-10000, -10000, -10000])

            subscore[0] = score1[0, j - 1, k - 1] + scoredArray[4, seqLen2[jStart + j], seqLen3[kStart + k]]
            subscore[1] = score1[0, j - 1, k] + scoredArray[4, seqLen2[jStart + j], 4]
            subscore[2] = score1[0, j, k - 1] + scoredArray[4, 4, seqLen3[kStart + k]]

            scoreXYplane = subscore[0]

            for i in range(1, 3):
                if subscore[i] > scoreXYplane:
                    scoreXYplane = subscore[i]

            score1[0, j, k] = scoreXYplane

    lengthy = 0

    for i in range(1, (mid - iStart + 1)):

        lengthy = (lengthy + 1) % 2
        score1[lengthy, 0, 0] = score1[(lengthy + 1) % 2, 0, 0] + scoredArray[seqLen1[iStart + i], 4, 4]

        for x in range(1, lengthj + 1):
            scoreXYplane = 0
            subscore2 = np.array([-10000, -10000, -10000])

            subscore2[0] = score1[(lengthy + 1)%2, x - 1, 0] + scoredArray[seqLen1[iStart + i], seqLen2[jStart + x], 4]
            subscore2[1] = score1[(lengthy + 1)%2, x, 0] + scoredArray[seqLen1[iStart + i], 4, 4]
            subscore2[2] = score1[lengthy, x - 1, 0] + scoredArray[4, seqLen2[jStart + x], 4]

            scoreXYplane = subscore2[0]

            for j in range(1, 3):
                if subscore2[j] > scoreXYplane:
                    scoreXYplane = subscore2[j]

            score1[lengthy, x, 0] = scoreXYplane

        for z in range(1, lengthk + 1):

            subscore3 = np.array([-10000, -10000, -10000])
            scoreXYplane = 0

            subscore3[0] = score1[((lengthy + 1) % 2), 0, z - 1] + scoredArray[seqLen1[iStart + i], 4, seqLen3[kStart + z]]
            subscore3[1] = score1[((lengthy + 1) % 2), 0, z] + scoredArray[seqLen1[iStart + i], 4, 4]
            subscore3[2] = score1[lengthy, 0, z - 1] + scoredArray[4, 4, seqLen3[kStart + z]]

            scoreXYplane = subscore3[0]

            for j in range(1, 3):
                if subscore3[j] > scoreXYplane:
                    scoreXYplane = subscore3[j]

            score1[lengthy, 0, z] = scoreXYplane

        for j in range(1, lengthj):
            for k in range(1, lengthk):
                subscore4 = np.array([-10000, -10000, -10000, -10000, -10000, -10000, -10000])
                scoreXYplane = 0

                subscore4[0] = score1[(lengthy + 1)%2, j - 1, k - 1] + scoredArray[seqLen1[iStart + i], seqLen2[jStart + j], seqLen3[kStart + k]]
                subscore4[1] = score1[(lengthy + 1)%2, j - 1, k] + scoredArray[seqLen1[iStart + i], seqLen2[jStart + j], 4]
                subscore4[2] = score1[(lengthy + 1)%2, j, k -1] + scoredArray[seqLen1[iStart + i], 4, seqLen3[kStart + k]]
                subscore4[3] = score1[lengthy, j -1, k - 1] + scoredArray[4, seqLen2[jStart +j], seqLen3[kStart + k]]
                subscore4[4] = score1[(lengthy + 1)%2, j, k] + scoredArray[seqLen1[iStart + i], 4, 4]
                subscore4[5] = score1[lengthy, j - 1, k] + scoredArray[4, seqLen2[jStart + j], 4]
                subscore4[6] = score1[lengthy, j, k - 1] + scoredArray[4, 4, seqLen3[kStart + k]]

                scoreXYplane = subscore4[0]

                for w in range(1, 7):
                    if subscore4[w] > scoreXYplane:
                        scoreXYplane = subscore4[w]

                score1[lengthy, j, k] = scoreXYplane

    forwardScore = lengthy

    #Score Backward

    score2[0, lengthj, lengthk] = 0
    for k in range(lengthk, 0, -1):
        score2[0, lengthj, k] = scoredArray[4, 4, seqLen3[kStart + k]] * (lengthk - k)

    for j in range(lengthj, 0,-1):
        score2[0, j, lengthk] = scoredArray[4, seqLen2[jStart + j], 4] * (lengthj - j)

    for j in range(lengthj - 1, 0, -1):
        for k in range(lengthk - 1, 0, -1):
            scoreXYplane = 0
            subscore5 = [-10000, -10000, -10000]

            subscore5[0] = score2[0, j + 1, k + 1] + scoredArray[4, seqLen2[jStart + j + 1], seqLen3[kStart + k + 1]]
            subscore5[1] = score2[0, j + 1, k] + scoredArray[4, seqLen2[jStart + j + 1], 4]
            subscore5[2] = score2[0, j, k + 1] + scoredArray[4, 4, seqLen3[kStart + k + 1]]

            scoreXYplane = subscore5[0]

            for i in range(1, 3):
                if subscore5[i] > scoreXYplane:
                    scoreXYplane = subscore5[i]

            score2[0, j, k] = scoreXYplane

    lengthy = 0

    for i in range(lengthi - 1, mid - iStart - 1, -1):

        lengthy = (lengthy + 1) % 2

        score2[lengthy, lengthj, lengthk] = score2[(lengthy + 1) % 2, lengthj, lengthk] + scoredArray[seqLen1[iStart + i + 1], 4, 4]
        for x in range(lengthj - 1, -1, -1):
            subscore6 = np.array([-10000, -10000, -10000])
            scoreXYplane = 0

            subscore6[0] = score2[(lengthy + 1) % 2, x + 1, lengthk] + scoredArray[seqLen1[iStart + i + 1], seqLen2[jStart + x + 1], 4]
            subscore6[1] = score2[(lengthy + 1) % 2, x, lengthk] + scoredArray[seqLen1[iStart + i + 1], 4, 4]
            subscore6[2] = score2[lengthy, x + 1, lengthk] + scoredArray[4, seqLen2[jStart + x + 1], 4]

            scoreXYplane = subscore6[0]

            for i in range(1, 3):
                if subscore6[i] > scoreXYplane:
                    scoreXYplane = subscore6[i]


            score2[lengthy, x, lengthk] = scoreXYplane

        for x in range(lengthk - 1, -1, -1):
            subscore7 = np.array([-10000, -10000, -10000])
            scoreXYplane = 0

            subscore7[0] = score2[(lengthy + 1) % 2, lengthj, x + 1] + scoredArray[seqLen1[iStart + i + 1], 4, seqLen3[kStart + x + 1]]
            subscore7[1] = score2[(lengthy + 1) % 2, lengthj, x] + scoredArray[seqLen1[iStart + i + 1], 4, 4]
            subscore7[2] = score2[lengthy, lengthj, x + 1] + scoredArray[4, 4, seqLen3[kStart + x + 1]]

            scoreXYplane = subscore7[0]

            for i in range(1, 3):
                if subscore7[i] > scoreXYplane:
                    scoreXYplane = subscore7[i]

            score2[lengthy, lengthj, x] = scoreXYplane

        for j in range(lengthj - 1, -1, -1):
            for k in range(lengthk - 1, -1, -1):

                subscore8 = np.array([-10000, -10000, -10000, -10000, -10000, -10000, -10000])

                scoreXYplane = 0

                subscore8[0] = score2[(lengthy + 1) % 2, j + 1, k + 1] + scoredArray[seqLen1[iStart + i + 1], seqLen2[jStart + j + 1], seqLen3[kStart + k + 1]]
                subscore8[1] = score2[(lengthy + 1) % 2, j + 1, k] + scoredArray[seqLen1[iStart + i + 1], seqLen2[jStart + j + 1], 4]
                subscore8[2] = score2[(lengthy + 1) % 2, j, k + 1] + scoredArray[seqLen1[iStart + i + 1], 4, seqLen3[kStart + k + 1]]
                subscore8[3] = score2[lengthy, j + 1, k + 1] + scoredArray[4, seqLen2[jStart + j + 1], seqLen3[kStart + k + 1]]
                subscore8[4] = score2[(lengthy + 1) % 2, j + 1, k] + scoredArray[seqLen1[iStart + i + 1], 4, 4]
                subscore8[5] = score2[lengthy, j + 1, k] + scoredArray[4, seqLen2[jStart + j + 1], 4]
                subscore8[6] = score2[lengthy, j, k + 1] + scoredArray[4, 4, seqLen3[kStart + k + 1]]

                scoreXYplane = subscore8[0]

                for i in range(1, 7):
                    if subscore8[i] > scoreXYplane:
                        scoreXYplane = subscore8[i]

                score2[lengthy, j, k] = scoreXYplane

    backwardScore = lengthy

    maxScore = score1[forwardScore, 0 , 0] + score2[backwardScore, 0, 0]

    print ("This is mid:" + str(mid))
    midX = mid
    midY = jStart
    midZ = kStart

    for j in range(0, lengthj):
        for k in range(0, lengthk):
            if (score1[forwardScore, j, k] + score2[backwardScore, j, k]) > maxScore:
                maxScore = score1[forwardScore, j, k] + score2[backwardScore, j, k]
                midY = jStart + j
                midZ = kStart + k

    global maxValue

    if maxValue:
        maxValue = False
        print ("Alignment Score:" + str(maxScore))

    alignment(iStart, jStart, kStart, midX, midY, midZ)
    alignment(midX, midY, midZ, LengthSeq1, LengthSeq2, LengthSeq3)

    return

def sequenceReturnNumbers(i):

    if str(i) == "A":
        return 0
    elif str(i) == "T":
        return 1
    elif str(i) == "G":
        return 2
    elif str(i) == "C":
        return 3
    else:
        return 4


def main():

    global MATCH
    MATCH = 5

    global MISMATCH
    MISMATCH = -4

    global INDEL
    INDEL = -8

    global scoredArray
    scoredArray = createScoring()

    global maxValue
    maxValue = True

    sequence1 = open("set1Seq1.txt", "r")
    sequence2 = open("set1Seq2.txt", "r")
    sequence3 = open("set1Seq3.txt", "r")

    sequence1Lines = sequence1.readlines()
    sequence2Lines = sequence2.readlines()
    sequence3Lines = sequence3.readlines()

    global seqLen1
    global seqLen2
    global seqLen3

    seqLen1Again = []
    seqLen2Again = []
    seqLen3Again = []

    for i in sequence1Lines:
        for letter in i:
            seqLen1Again.append(sequenceReturnNumbers(letter))

    for i in sequence2Lines:
        for letter in i:
            seqLen2Again.append(sequenceReturnNumbers(letter))

    for i in sequence3Lines:
        for letter in i:
            seqLen3Again.append(sequenceReturnNumbers(letter))

    seqLen1 = np.array(seqLen1Again)
    seqLen2 = np.array(seqLen2Again)
    seqLen3 = np.array(seqLen3Again)

    alignment(0, 0, 0, len(seqLen1) - 1, len(seqLen2) - 1, len(seqLen3) - 1)

main()
