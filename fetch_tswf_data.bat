@echo off

echo %time%

call activate nanodust

(
	start python solo_download_wrap.py "20200615" "20201231"
	start python solo_download_wrap.py "20210101" "20210630"
	start python solo_download_wrap.py "20210701" "20211231"
	start python solo_download_wrap.py "20220101" "20220630"
	start python solo_download_wrap.py "20220701" "20201231"
	start python solo_download_wrap.py "20230101" "20230630"
	
) | set /P "="

echo %time%
