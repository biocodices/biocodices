set :application, 'biocodices'
set :repo_url, 'git@github.com:biocodices/biocodices.git'
set :deploy_to, '~/biocodices/app'
set :keep_releases, 3

namespace :deploy do
  task :install_app do
    on roles(:app) do
      conda_path = "~/anaconda3/bin"
      python = File.join(conda_path, "python")
      execute "cd #{current_path} && #{python} setup.py install"
      execute "export PATH=$PATH:#{conda_path} && bioco --version"
    end
  end

  after :deploy, "deploy:install_app"
end
