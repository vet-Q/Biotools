package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

func Readfasta(filePath string) (result []string) {

	file, err := os.Open(filePath)
	if err != nil {

		// 파일이 존재하지 않는 경우 에러를 반환하는 코드
		if os.IsNotExist(err) {
			fmt.Printf("파일이 존재하지 않습니다: %v", err)
		} else {
			fmt.Printf("파일을 열 수 없습니다: %v", err)
		}
		return
	}
	defer file.Close()

	scanner := bufio.NewReader(file)

	for {
		line, err := scanner.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			fmt.Printf("파일을 읽는 도중 오류가 발생했습니다: %v", err)
			return
		}

		// 각 줄을 처리하는 로직을 구현합니다.
		if strings.HasPrefix(line, ">") {
			result = append(result, line)
		}

	}
	return result
}

func Writetabular(records []string, inputfilePath string, outputfilePath string) error {
	file, err := os.Open(inputfilePath)
	if err != nil {

		// 파일이 존재하지 않는 경우 에러를 반환하는 코드
		if os.IsNotExist(err) {
			fmt.Printf("파일이 존재하지 않습니다: %v", err)
		} else {
			fmt.Printf("파일을 열 수 없습니다: %v", err)
		}
		return nil
	}
	defer file.Close()

	scanner := bufio.NewReader(file)

	for {
		line, err := scanner.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			fmt.Printf("파일을 읽는 도중 오류가 발생했습니다: %v", err)
			return nil
		}

		// 각 줄을 처리하는 로직을 구현합니다.
		if strings.HasPrefix(line, ">") {
			result = append(result, line)
		}

	}

	outfile, err := os.Create(outputfilePath)
	if err != nil {
		return err
	}
	defer outfile.Close()

	writer := bufio.NewWriter(outfile)
	for _, record := range records {
		_, err := fmt.Fprintf(writer, ">%s\n%s\n", record.Header, record.Seq)
		if err != nil {
			return err
		}
	}

	err = writer.Flush()
	if err != nil {
		return err
	}

	return nil
}

func main() {
	filePath := "C:\\Users\\user\\Desktop\\ASF_KNU\\upload_ASF\\Galaxy2-[Genbank_to_Five_Column_Format_on_data_1] (1).tabular"
	// readfilePath := "C:\\Users\\user\\Desktop\\ASF_KNU\\ASF_knu_final.fasta.fasta"
	fmt.Println(Readfasta(filePath))
}
