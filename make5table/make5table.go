package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strings"
)

func main() {
	// 원본 FASTA 파일 경로
	originalFilePath := "C:\\Users\\kwono\\Desktop\\ASF_feature_table_all.tabular"

	// 복사 대상 FASTA 파일 경로
	targetFilePath := "C:\\Users\\kwono\\Desktop\\ASF_feature_table_all.tabular_modified"

	slicedPath := "C:\\Users\\kwono\\Desktop\\merged.fasta"
	// 슬라이스된 이름들
	names := slicedData(slicedPath)

	err := copyFastaFile(originalFilePath, targetFilePath, names)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("FASTA 파일이 성공적으로 복사되었습니다.")
}

// siice된 이름을 추출하는 함수
func slicedData(inputData string) (result []string) {
	originalFile, err := os.Open(inputData)
	if err != nil {
		return
	}
	defer originalFile.Close()

	// 파일 스캐너 생성
	scanner := bufio.NewScanner(originalFile)

	for scanner.Scan() {
		line := scanner.Text()

		// 헤더인 경우 이름 변경
		if strings.HasPrefix(line, ">") {
			line = strings.TrimPrefix(line, ">")
			result = append(result, line)
		}
	}
	return result
}

// FASTA 파일 복사 함수
func copyFastaFile(originalFilePath, targetFilePath string, names []string) error {
	// 원본 파일 열기
	originalFile, err := os.Open(originalFilePath)
	if err != nil {
		return err
	}
	defer originalFile.Close()

	// 복사 대상 파일 생성
	targetFile, err := os.Create(targetFilePath)
	if err != nil {
		return err
	}
	defer targetFile.Close()

	// 파일 스캐너 생성
	scanner := bufio.NewScanner(originalFile)

	// 대상 파일 작성기 생성
	writer := bufio.NewWriter(targetFile)

	// FASTA 파일 복사
	err = processFastaData(scanner, writer, names)
	if err != nil {
		return err
	}

	// 대상 파일 버퍼 비우고 저장
	err = writer.Flush()
	if err != nil {
		return err
	}

	return nil
}

// FASTA 데이터 처리 함수
func processFastaData(scanner *bufio.Scanner, writer *bufio.Writer, names []string) error {
	for scanner.Scan() {
		line := scanner.Text()

		// 헤더인 경우 이름 변경
		if strings.HasPrefix(line, ">") {
			header := strings.TrimPrefix(line, ">")
			newHeader := getNewHeader(header, names)

			// 대상 파일에 쓰기
			_, err := writer.WriteString(">" + newHeader + "\n")
			if err != nil {
				return err
			}
		} else {
			// 데이터인 경우 그대로 대상 파일에 쓰기
			_, err := writer.WriteString(line + "\n")
			if err != nil {
				return err
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return err
	}

	return nil
}

// 슬라이스된 이름을 가져오는 함수
func getNewHeader(oldHeader string, names []string) string {
	// 예외 처리: 이름이 없는 경우 기본값 사용
	if len(names) == 0 {
		return "NewName"
	}

	// 이름 슬라이스 인덱스 가져오기
	index := len(names) - 1

	// 슬라이스된 이름 반환
	return names[index]
}
